//! Queries a BAM file with a given region.
//!
//! The input BAM must have an index in the same directory.
//!
//! The result matches the output of `samtools view <src> <region>`.

use std::{env, fs::File, io};

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::io::Write};

const UNMAPPED: &str = "*";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let raw_region = args.next().expect("missing region");

    let mut reader = File::open(&src).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;

    let index = bam::fs::read_associated_index(&src)?;

    let records: Box<dyn Iterator<Item = io::Result<bam::Record>>> = if raw_region == UNMAPPED {
        reader.query_unmapped(&index).map(Box::new)?
    } else {
        let region = raw_region.parse()?;
        Box::new(reader.query(&header, &index, &region)?.records())
    };

    let stdout = io::stdout().lock();
    let mut writer = sam::io::Writer::new(stdout);

    for result in records {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    writer.finish(&header)?;

    Ok(())
}
