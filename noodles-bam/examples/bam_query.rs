//! Querys a BAM file with a given region.
//!
//! The input BAM must have an index in the same directory.
//!
//! While the results are not formatted the same, the records printed match the output of `samtools
//! view <src> <region>`.

use std::{env, fs::File, path::PathBuf, str};

use noodles_bam::{self as bam, bai};
use noodles_sam as sam;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args();

    let src = args.nth(1).map(PathBuf::from).expect("missing src");
    let region = args.next().expect("missing region").parse()?;

    let mut reader = File::open(&src).map(bam::Reader::new)?;
    let header: sam::Header = reader.read_header()?.parse()?;
    let reference_sequences = header.reference_sequences();

    let index = bai::read(src.with_extension("bam.bai"))?;

    let query = reader.query(reference_sequences, &index, &region)?;

    for result in query {
        let record = result?;

        let name = str::from_utf8(record.name())?;

        let reference_sequence_id = record.reference_sequence_id() as usize;
        let (_, reference_sequence) = reference_sequences
            .get_index(reference_sequence_id)
            .ok_or_else(|| "invalid reference sequence id")?;

        let start = record.position() + 1;
        let len = record.cigar().mapped_len() as i32;
        let end = start + len - 1;

        println!("{} ({}:{}-{})", name, reference_sequence.name(), start, end);
    }

    Ok(())
}
