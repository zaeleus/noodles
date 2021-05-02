//! Prints a BAM file in the SAM format.
//!
//! The result matches the output of `samtools view <src>`.

use std::{env, fs::File};

use noodles_bam as bam;
use noodles_sam as sam;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bam::Reader::new)?;
    let header: sam::Header = reader.read_header()?.parse()?;
    reader.read_reference_sequences()?;

    for result in reader.records() {
        let record = result?;
        let sam_record = record.try_into_sam_record(header.reference_sequences())?;
        println!("{}", sam_record);
    }

    Ok(())
}
