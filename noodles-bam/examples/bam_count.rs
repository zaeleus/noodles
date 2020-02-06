//! Counts the number of records in a BAM file.
//!
//! The result matches the output of `samtools view -c <src>`.

use std::{env, fs::File};

use noodles_bam as bam;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let file = File::open(src)?;
    let mut reader = bam::Reader::new(file);
    reader.header()?;

    let count = reader.records().count();
    println!("{}", count);

    Ok(())
}
