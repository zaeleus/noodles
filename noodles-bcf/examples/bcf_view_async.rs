//! Prints a BCF file in the VCF format.
//!
//! The result matches the output of `bcftools view <src>`.

use std::env;

use futures::TryStreamExt;
use noodles_bcf as bcf;
use tokio::fs::File;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(bcf::AsyncReader::new)?;
    reader.read_file_format().await?;

    let raw_header = reader.read_header().await?;
    let header = raw_header.parse()?;
    let string_maps = raw_header.parse()?;

    print!("{raw_header}");

    let mut records = reader.lazy_records();

    while let Some(record) = records.try_next().await? {
        let vcf_record = record.try_into_vcf_record(&header, &string_maps)?;
        println!("{vcf_record}");
    }

    Ok(())
}
