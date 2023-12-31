//! Queries a BAM file with a given region.
//!
//! The input BAM must have an index in the same directory.
//!
//! The result matches the output of `samtools view <src> <region>`.

use std::{env, io, path::PathBuf};

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::RecordBuf};
use sam::AlignmentWriter;

const UNMAPPED: &str = "*";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args();

    let src = args.nth(1).map(PathBuf::from).expect("missing src");
    let raw_region = args.next().expect("missing region");

    let mut reader = bam::indexed_reader::Builder::default().build_from_path(src)?;
    let header = reader.read_header()?;

    let records: Box<dyn Iterator<Item = io::Result<RecordBuf>>> = if raw_region == UNMAPPED {
        reader.query_unmapped(&header).map(Box::new)?
    } else {
        let region = raw_region.parse()?;
        reader.query(&header, &region).map(Box::new)?
    };

    let stdout = io::stdout().lock();
    let mut writer = sam::Writer::new(stdout);

    for result in records {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    writer.finish(&header)?;

    Ok(())
}
