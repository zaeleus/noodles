//! Counts the number of records in a BAM file.
//!
//! The result matches the output of `samtools view --count <src>`.

use std::env;

use futures::TryStreamExt;
use noodles_bam as bam;
use tokio::fs::File;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(bam::r#async::io::Reader::new)?;
    let header = reader.read_header().await?;

    let mut records = reader.record_bufs(&header);
    let mut n = 0;

    while records.try_next().await?.is_some() {
        n += 1;
    }

    println!("{n}");

    Ok(())
}
