//! Queries a bgzipped SAM file with a given region.
//!
//! The input bgzipped SAM file must have an associated coordinate-sorted index (CSI) in the same
//! directory.
//!
//! The result matches the output of `samtools view <src> <region>`.

use std::{env, io, path::PathBuf};

use noodles_sam as sam;

const UNMAPPED: &str = "*";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args();

    let src = args.nth(1).map(PathBuf::from).expect("missing src");
    let raw_region = args.next().expect("missing region");

    let mut reader = sam::io::indexed_reader::Builder::default().build_from_path(src)?;
    let header = reader.read_header()?;

    let records: Box<dyn Iterator<Item = io::Result<sam::Record>>> = if raw_region == UNMAPPED {
        reader.query_unmapped().map(Box::new)?
    } else {
        let region = raw_region.parse()?;
        reader
            .query(&header, &region)
            .map(|query| Box::new(query.records()))?
    };

    let stdout = io::stdout().lock();
    let mut writer = sam::io::Writer::new(stdout);

    for result in records {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    Ok(())
}
