//! Prints a BAM file in the SAM format.
//!
//! The result matches the output of `samtools view <src>`.

use std::env;

use futures::TryStreamExt;
use noodles_bam as bam;
use noodles_sam as sam;
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(bam::AsyncReader::new)?;
    let header: sam::Header = reader.read_header().await?.parse()?;
    reader.read_reference_sequences().await?;

    let mut records = reader.records();

    let mut writer = sam::AsyncWriter::new(io::stdout());

    while let Some(record) = records.try_next().await? {
        writer.write_alignment_record(&header, &record).await?;
    }

    Ok(())
}
