//! Creates a new VCF file.
//!
//! This writes a VCF header and one (sparse) record to stdout.
//!
//! Verify the output by piping to `bcftools view --no-version`.

use noodles_vcf::{self as vcf, header::Contig, record::Position};
use tokio::io;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = vcf::AsyncWriter::new(io::stdout());

    let header = vcf::Header::builder()
        .add_contig(Contig::new(String::from("sq0")))
        .build();

    writer.write_header(&header).await?;

    let record = vcf::Record::builder()
        .set_chromosome("sq0".parse()?)
        .set_position(Position::try_from(1)?)
        .set_reference_bases("A".parse()?)
        .build()?;

    writer.write_record(&record).await?;

    Ok(())
}
