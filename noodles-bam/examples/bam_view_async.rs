//! Prints a BAM file in the SAM format.
//!
//! The result matches the output of `samtools view <src>`.

use std::env;

use futures::TryStreamExt;
use noodles_bam as bam;
use noodles_sam as sam;
use tokio::{
    fs::File,
    io::{self, AsyncWriteExt},
};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(bam::r#async::io::Reader::new)?;
    let header = reader.read_header().await?;

    let mut records = reader.record_bufs(&header);

    let mut writer = sam::r#async::io::Writer::new(io::stdout());

    while let Some(record) = records.try_next().await? {
        writer.write_alignment_record(&header, &record).await?;
    }

    writer.get_mut().shutdown().await?;

    Ok(())
}
