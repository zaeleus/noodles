//! Validates and prints a VCF file to stdout.
//!
//! The result matches the output of `bcftools view <src>`.

use std::env;

use futures::TryStreamExt;
use noodles_vcf as vcf;
use tokio::{
    fs::File,
    io::{self, BufReader},
};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .await
        .map(BufReader::new)
        .map(vcf::AsyncReader::new)?;

    let header = reader.read_header().await?;

    let mut writer = vcf::AsyncWriter::new(io::stdout());
    writer.write_header(&header).await?;

    let mut records = reader.records(&header);

    while let Some(record) = records.try_next().await? {
        writer.write_record(&record).await?;
    }

    Ok(())
}
