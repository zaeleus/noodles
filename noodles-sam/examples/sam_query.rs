//! Queries a bgzipped SAM file with a given region.
//!
//! The input bgzipped SAM file must have an associated coordinate-sorted index (CSI) in the same
//! directory.
//!
//! The result matches the output of `samtools view <src> <region>`.

use std::{env, fs::File, io};

use noodles_bgzf as bgzf;
use noodles_sam as sam;

const UNMAPPED: &str = "*";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let raw_region = args.next().expect("missing region");

    let mut reader = File::open(&src)
        .map(bgzf::io::Reader::new)
        .map(sam::io::Reader::new)?;

    let header = reader.read_header()?;

    let index = sam::fs::read_associated_index(&src)?;

    let records: Box<dyn Iterator<Item = io::Result<sam::Record>>> = if raw_region == UNMAPPED {
        reader.query_unmapped(&index).map(Box::new)?
    } else {
        let region = raw_region.parse()?;

        reader
            .query(&header, &index, &region)
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
