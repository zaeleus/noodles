//! Prints the header of a BAM file.
//!
//! A BAM file header is a SAM header.
//!
//! The result matches the output of `samtools view --no-PG -H <src>`.

use std::{env, fs::File, io};

use noodles_bam as bam;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");
    let mut reader = File::open(src).map(bam::Reader::new)?;

    let header = reader.read_header()?;
    print!("{}", header);

    Ok(())
}
