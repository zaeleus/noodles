//! Queries a variant file with the given region.
//!
//! The input must have an associated index in the same directory.
//!
//! The result matches the output of `bcftools view <src> <region>`.

use std::env;

use futures::TryStreamExt;
use noodles_util::variant;
use noodles_vcf as vcf;
use tokio::io::{self, AsyncWriteExt};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let region = args.next().expect("missing region").parse()?;

    let mut reader = variant::r#async::io::reader::Builder::default()
        .build_from_path(&src)
        .await?;

    let header = reader.read_header().await?;

    let index = variant::r#async::fs::read_associated_index(&src).await?;
    let mut query = reader.query(&header, &index, &region)?;

    let mut writer = vcf::r#async::io::Writer::new(io::stdout());

    writer.write_header(&header).await?;

    while let Some(record) = query.try_next().await? {
        writer
            .write_variant_record(&header, record.as_ref())
            .await?;
    }

    writer.get_mut().shutdown().await?;

    Ok(())
}
