//! Prints a BAM file in the SAM format using multiple threads.
//!
//! The result matches the output of `samtools view <src>`.

use std::{env, fs::File, io};

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::io::Write};
use noodles_bgzf as bgzf;

fn main() -> io::Result<()> {
    let reader_threads = 2;
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(|file| bgzf::MultithreadedReader::with_worker_count(std::num::NonZeroUsize::new(reader_threads).unwrap(), file))
        .map(bam::io::Reader::from)?;

    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = sam::io::Writer::new(stdout);

    for result in reader.records() {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    writer.finish(&header)?;

    Ok(())
}
