//! Queries a BAM file with a given region.
//!
//! The input BAM must have an index in the same directory.
//!
//! The result matches the output of `samtools view <src> <region>`.

use std::{env, fs::File, path::PathBuf};

use noodles_bam::{self as bam, bai};
use noodles_core::Region;
use noodles_sam as sam;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args();

    let src = args.nth(1).map(PathBuf::from).expect("missing src");
    let raw_region = args.next().expect("missing region");

    let mut reader = File::open(&src).map(bam::Reader::new)?;
    let header: sam::Header = reader.read_header()?.parse()?;
    let reference_sequences = header.reference_sequences();

    let index = bai::read(src.with_extension("bam.bai"))?;

    let region = Region::from_str_reference_sequences(&raw_region, reference_sequences)?;
    let query = reader.query(reference_sequences, &index, &region)?;

    for result in query {
        let record = result?;
        let sam_record = record.try_into_sam_record(reference_sequences)?;
        println!("{}", sam_record);
    }

    Ok(())
}
