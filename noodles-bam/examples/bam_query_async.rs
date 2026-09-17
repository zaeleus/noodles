//! Queries a BAM file with a given region.
//!
//! The input BAM must have an index in the same directory.
//!
//! The result matches the output of `samtools view <src> <region>`.

use std::{env, pin::Pin};

use futures::{Stream, TryStreamExt};
use noodles_bam as bam;
use noodles_sam as sam;
use tokio::{fs::File, io};

const UNMAPPED: &str = "*";

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let raw_region = args.next().expect("missing region");

    let mut reader = File::open(&src).await.map(bam::r#async::io::Reader::new)?;
    let header = reader.read_header().await?;

    let index = bam::r#async::fs::read_associated_index(&src).await?;

    let mut records: Pin<Box<dyn Stream<Item = io::Result<bam::Record>>>> =
        if raw_region == UNMAPPED {
            reader.query_unmapped(&index).await.map(Box::pin)?
        } else {
            let region = raw_region.parse()?;
            Box::pin(reader.query(&header, &index, &region)?.records())
        };

    let mut writer = sam::r#async::io::Writer::new(io::stdout());

    while let Some(record) = records.try_next().await? {
        writer.write_alignment_record(&header, &record).await?;
    }

    Ok(())
}
