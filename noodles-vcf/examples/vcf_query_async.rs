//! Queries a VCF file with a given region.
//!
//! The input VCF must have an associated index in the same directory.
//!
//! The result matches the output `bcftools view --no-header <src> <region>`.

use std::env;

use futures::TryStreamExt;
use noodles_bgzf as bgzf;
use noodles_vcf as vcf;
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let region = args.next().map(|s| s.parse()).expect("missing region")?;

    let mut reader = File::open(&src)
        .await
        .map(bgzf::r#async::io::Reader::new)
        .map(vcf::r#async::io::Reader::new)?;

    let header = reader.read_header().await?;

    let index = vcf::r#async::fs::read_associated_index(&src).await?;
    let mut query = reader.query(&header, &index, &region)?.records();

    let mut writer = vcf::r#async::io::Writer::new(io::stdout());

    while let Some(record) = query.try_next().await? {
        writer.write_variant_record(&header, &record).await?;
    }

    writer.shutdown().await?;

    Ok(())
}
