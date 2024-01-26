//! Prints the header of a BCF file.
//!
//! The result matches the output of `bcftools head <src>`.

use std::env;

use noodles_bcf as bcf;
use noodles_vcf as vcf;
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(bcf::r#async::io::Reader::new)?;
    reader.read_file_format().await?;
    let header = reader.read_header().await?;

    let mut writer = vcf::r#async::io::Writer::new(io::stdout());
    writer.write_header(&header).await?;

    writer.shutdown().await?;

    Ok(())
}
