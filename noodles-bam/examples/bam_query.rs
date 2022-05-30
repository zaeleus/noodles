//! Queries a BAM file with a given region.
//!
//! The input BAM must have an index in the same directory.
//!
//! The result matches the output of `samtools view <src> <region>`.

use std::{env, fs::File, io, path::PathBuf};

use noodles_bam::{self as bam, bai};
use noodles_sam as sam;
use sam::AlignmentWriter;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args();

    let src = args.nth(1).map(PathBuf::from).expect("missing src");
    let region = args.next().expect("missing region").parse()?;

    let mut reader = File::open(&src).map(bam::Reader::new)?;

    let header: sam::Header = reader.read_header()?.parse()?;
    let index = bai::read(src.with_extension("bam.bai"))?;
    let query = reader.query(header.reference_sequences(), &index, &region)?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = sam::Writer::new(handle);

    for result in query {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    writer.finish(&header)?;

    Ok(())
}
