//! Prints a variant file in the VCF format.
//!
//! The result matches the output of `bcftools view <src>`.

use std::env;

use futures::TryStreamExt;
use noodles_util::variant;
use noodles_vcf as vcf;
use tokio::{
    fs::File,
    io::{self, AsyncRead, AsyncWriteExt},
};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let source: Box<dyn AsyncRead + Unpin> = if src == "-" {
        Box::new(io::stdin())
    } else {
        File::open(src).await.map(Box::new)?
    };

    let mut reader = variant::r#async::io::Reader::new(source).await?;
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
