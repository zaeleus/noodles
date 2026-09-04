//! Converts a SAM file to BAM, encoding records on a thread pool.
//!
//! `bgzf::io::MultithreadedWriter` compresses on several threads, but the records themselves are
//! parsed and encoded by whichever thread calls the writer. Past the point where compression keeps
//! up, that thread is the constraint, and adding compression workers does nothing.
//!
//! This reads the input in chunks, splits each chunk into one contiguous range per worker aligned
//! to record boundaries, encodes the ranges in parallel, and writes the results in order. Only the
//! boundary rule is specific to SAM: records are one line, so a range ends at a line feed.
//!
//! The output is identical to writing the records one at a time.
//!
//! The equivalent of `samtools view -b -o <dst> <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufReader, Read, Write},
    ops::Range,
};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::{self as sam, alignment::io::Write as _};
use rayon::prelude::*;

const CHUNK_SIZE: usize = 32 << 20;
const LINE_FEED: u8 = b'\n';

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");
    let dst = args.next().expect("missing dst");

    let worker_count = std::thread::available_parallelism()?;

    let mut reader = sam::io::Reader::new(BufReader::new(File::open(src)?));
    let header = reader.read_header()?;

    let encoder =
        bgzf::io::MultithreadedWriter::with_worker_count(worker_count, File::create(dst)?);
    let mut writer = bam::io::Writer::from(encoder);
    writer.write_header(&header)?;

    let mut inner = reader.into_inner();
    let mut buf = vec![0; CHUNK_SIZE];
    let mut carry_len = 0;

    loop {
        let mut filled = carry_len;

        while filled < buf.len() {
            let n = inner.read(&mut buf[filled..])?;

            if n == 0 {
                break;
            }

            filled += n;
        }

        if filled == 0 {
            break;
        }

        // Only complete records are handed out. The tail moves to the front for the next round.
        let end = match buf[..filled].iter().rposition(|&b| b == LINE_FEED) {
            Some(i) => i + 1,
            None => filled,
        };

        let ranges = split(&buf[..end], worker_count.get());

        let blocks: Vec<Vec<u8>> = ranges
            .par_iter()
            .map(|range| encode(&header, &buf[range.clone()]))
            .collect::<io::Result<_>>()?;

        for block in &blocks {
            writer.get_mut().write_all(block)?;
        }

        buf.copy_within(end..filled, 0);
        carry_len = filled - end;

        if filled < CHUNK_SIZE {
            break;
        }
    }

    writer.finish(&header)?;

    Ok(())
}

/// Splits a buffer of complete records into `n` contiguous ranges, each ending on a record
/// boundary.
fn split(src: &[u8], n: usize) -> Vec<Range<usize>> {
    let mut ranges = Vec::with_capacity(n);
    let target = (src.len() / n).max(1);
    let mut start = 0;

    while start < src.len() {
        let mut end = (start + target).min(src.len());

        if end < src.len() {
            match src[end..].iter().position(|&b| b == LINE_FEED) {
                Some(i) => end += i + 1,
                None => end = src.len(),
            }
        }

        ranges.push(start..end);
        start = end;
    }

    ranges
}

fn encode(header: &sam::Header, src: &[u8]) -> io::Result<Vec<u8>> {
    let mut writer = bam::io::Writer::from(Vec::with_capacity(src.len()));

    for line in src.split(|&b| b == LINE_FEED) {
        if line.is_empty() {
            continue;
        }

        let record = sam::Record::try_from(line)?;
        writer.write_alignment_record(header, &record)?;
    }

    Ok(writer.into_inner())
}
