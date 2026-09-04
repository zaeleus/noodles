//! Converts a BAM file to SAM, formatting records on a thread pool.
//!
//! This is the counterpart of `bam_write_parallel`. `bgzf::io::MultithreadedReader` decompresses
//! on several threads, but the records are formatted into text by whichever thread drives the
//! writer, and past a few workers that thread is the constraint.
//!
//! Records are read into a batch, the batch is formatted in parallel, and the results are written
//! in order. Unlike the SAM input case there are no byte offsets to find: a batch of records is
//! already a list of independent units, which is what makes this the simpler direction.
//!
//! The output is identical to writing the records one at a time.
//!
//! The equivalent of `samtools view -h -o <dst> <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufWriter, Write},
};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::{self as sam, alignment::io::Write as _};
use rayon::prelude::*;

const BATCH_SIZE: usize = 1 << 16;

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");
    let dst = args.next().expect("missing dst");

    let worker_count = std::thread::available_parallelism()?;

    let decoder = bgzf::io::MultithreadedReader::with_worker_count(worker_count, File::open(src)?);
    let mut reader = bam::io::Reader::from(decoder);
    let header = reader.read_header()?;

    let mut inner = BufWriter::new(File::create(dst)?);

    {
        let mut writer = sam::io::Writer::new(&mut inner);
        writer.write_alignment_header(&header)?;
    }

    let mut batch = vec![bam::Record::default(); BATCH_SIZE];

    loop {
        let mut filled = 0;

        while filled < BATCH_SIZE {
            if reader.read_record(&mut batch[filled])? == 0 {
                break;
            }

            filled += 1;
        }

        if filled == 0 {
            break;
        }

        let chunk_size = (filled / worker_count).max(1);

        let blocks: Vec<Vec<u8>> = batch[..filled]
            .par_chunks(chunk_size)
            .map(|chunk| format(&header, chunk))
            .collect::<io::Result<_>>()?;

        for block in &blocks {
            inner.write_all(block)?;
        }

        if filled < BATCH_SIZE {
            break;
        }
    }

    Ok(())
}

fn format(header: &sam::Header, records: &[bam::Record]) -> io::Result<Vec<u8>> {
    let mut writer = sam::io::Writer::new(Vec::new());

    for record in records {
        writer.write_alignment_record(header, record)?;
    }

    Ok(writer.into_inner())
}
