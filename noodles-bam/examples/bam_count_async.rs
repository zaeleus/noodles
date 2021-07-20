//! Counts the number of records in a BAM file.
//!
//! The result matches the output of `samtools view -c <src>`.

use std::env;

use futures::StreamExt;
use noodles_bam as bam;
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(bam::AsyncReader::new)?;
    reader.read_header().await?;
    reader.read_reference_sequences().await?;

    let mut records = reader.records();
    let mut n = 0;

    while let Some(result) = records.next().await {
        let _ = result?;
        n += 1;
    }

    println!("{}", n);

    Ok(())
}
