//! Prints a variant file in the VCF format.
//!
//! The result matches the output of `bcftools view <src>`.

use std::env;

use futures::TryStreamExt;
use noodles_util::variant;
use noodles_vcf as vcf;
use tokio::io::{self, AsyncWriteExt};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let builder = variant::r#async::io::reader::Builder::default();

    let mut reader = if src == "-" {
        builder.build_from_reader(io::stdin()).await?
    } else {
        builder.build_from_path(src).await?
    };

    let header = reader.read_header().await?;

    let mut writer = vcf::r#async::io::Writer::new(io::stdout());
    writer.write_header(&header).await?;

    let mut records = reader.records();

    while let Some(record) = records.try_next().await? {
        writer
            .write_variant_record(&header, record.as_ref())
            .await?;
    }

    writer.get_mut().shutdown().await?;

    Ok(())
}
