//! Counts the number of records in a VCF file.
//!
//! The result matches the output of `bcftools view --no-header <src> | wc -l`.

use std::env;

use futures::TryStreamExt;
use noodles_vcf as vcf;
use tokio::{fs::File, io::BufReader};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .await
        .map(BufReader::new)
        .map(vcf::r#async::io::Reader::new)?;

    let header = reader.read_header().await?;

    let mut records = reader.record_bufs(&header);
    let mut n = 0;

    while records.try_next().await?.is_some() {
        n += 1;
    }

    println!("{n}");

    Ok(())
}
