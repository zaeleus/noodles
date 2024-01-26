//! Queries a BCF file with a given region.
//!
//! The input BCF must have an index in the same directory.
//!
//! The result matches the output of `bcftools view --no-header <src> <region>`.

use std::{env, path::PathBuf};

use futures::TryStreamExt;
use noodles_bcf::{self as bcf, header::StringMaps};
use noodles_csi as csi;
use noodles_vcf as vcf;
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args();

    let src = args.nth(1).map(PathBuf::from).expect("missing src");
    let raw_region = args.next().expect("missing region");

    let mut reader = File::open(&src).await.map(bcf::r#async::io::Reader::new)?;
    reader.read_file_format().await?;
    let raw_header = reader.read_header().await?;

    let header: vcf::Header = raw_header.parse()?;
    let string_maps: StringMaps = raw_header.parse()?;

    let index = csi::r#async::read(src.with_extension("bcf.csi")).await?;

    let region = raw_region.parse()?;
    let mut query = reader.query(string_maps.contigs(), &index, &region)?;

    let mut writer = vcf::r#async::io::Writer::new(io::stdout());

    while let Some(record) = query.try_next().await? {
        let vcf_record = record.try_into_vcf_record(&header, &string_maps)?;
        writer.write_record(&vcf_record).await?;
    }

    writer.shutdown().await?;

    Ok(())
}
