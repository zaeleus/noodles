//! Queries a VCF file with a given region.
//!
//! The input VCF must have an associated index in the same directory.
//!
//! The result matches the output `bcftools view --no-header <src> <region>`.

use std::{env, path::PathBuf};

use futures::TryStreamExt;
use noodles_bgzf as bgzf;
use noodles_tabix as tabix;
use noodles_vcf as vcf;
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args();

    let src = args.nth(1).map(PathBuf::from).expect("missing src");
    let region = args.next().map(|s| s.parse()).expect("missing region")?;

    let mut reader = File::open(&src)
        .await
        .map(bgzf::AsyncReader::new)
        .map(vcf::r#async::io::Reader::new)?;

    let header = reader.read_header().await?;

    let index = tabix::r#async::read(src.with_extension("gz.tbi")).await?;
    let mut query = reader.query(&header, &index, &region)?;

    let mut writer = vcf::r#async::io::Writer::new(io::stdout());

    while let Some(record) = query.try_next().await? {
        writer.write_record(&record).await?;
    }

    writer.shutdown().await?;

    Ok(())
}
