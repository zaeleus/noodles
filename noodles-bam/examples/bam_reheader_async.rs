//! Replaces the SAM header of a BAM file.
//!
//! This is similar to the functionality of `samtools reheader`.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::env;

use futures::TryStreamExt;
use noodles_bam as bam;
use noodles_sam as sam;
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(bam::AsyncReader::new)?;
    let mut header: sam::Header = reader.read_header().await?.parse()?;
    reader.read_reference_sequences().await?;

    header.add_comment("a comment added by noodles-bam");

    let mut writer = bam::AsyncWriter::new(io::stdout());
    writer.write_header(&header).await?;
    writer
        .write_reference_sequences(header.reference_sequences())
        .await?;

    let mut records = reader.records();

    while let Some(record) = records.try_next().await? {
        writer.write_record(&header, &record).await?;
    }

    writer.shutdown().await?;

    Ok(())
}
