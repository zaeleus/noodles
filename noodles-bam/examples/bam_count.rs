//! Counts the number of records in a BAM file.
//!
//! The result matches the output of `samtools view -c <src>`.

use std::{env, fs::File, io};

use noodles_bam as bam;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bam::Reader::new)?;
    reader.read_header()?;
    reader.read_reference_sequences()?;

    let count = reader.records().count();
    println!("{}", count);

    Ok(())
}
