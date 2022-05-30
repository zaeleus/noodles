//! Queries a BAM file with a given region.
//!
//! The input BAM must have an index in the same directory.
//!
//! The result matches the output of `samtools view <src> <region>`.

use std::{env, path::PathBuf};

use futures::TryStreamExt;
use noodles_bam::{self as bam, bai};
use noodles_sam as sam;
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args();

    let src = args.nth(1).map(PathBuf::from).expect("missing src");
    let region = args.next().expect("missing region").parse()?;

    let mut reader = File::open(&src).await.map(bam::AsyncReader::new)?;

    let header: sam::Header = reader.read_header().await?.parse()?;
    let index = bai::r#async::read(src.with_extension("bam.bai")).await?;
    let mut query = reader.query(header.reference_sequences(), &index, &region)?;

    let mut writer = sam::AsyncWriter::new(io::stdout());

    while let Some(record) = query.try_next().await? {
        writer.write_alignment_record(&header, &record).await?;
    }

    Ok(())
}
