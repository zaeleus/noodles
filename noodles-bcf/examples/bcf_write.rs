//! Creates a new BCF file.
//!
//! This writes a BCF file format, VCF header, and a single VCF to stdout.
//!
//! Verify the output by piping to `bcftools view --no-version`.

use std::{convert::TryFrom, io};

use noodles_bcf as bcf;
use noodles_vcf::{self as vcf, header::Contig, record::Position};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let stdout = io::stdout();
    let handle = stdout.lock();

    let mut writer = bcf::Writer::new(handle);
    writer.write_file_format()?;

    let header = vcf::Header::builder()
        .add_filter(vcf::header::Filter::pass())
        .add_contig(Contig::new("sq0"))
        .build();

    writer.write_header(&header)?;

    let raw_header = header.to_string();
    let string_map = raw_header.parse()?;

    let record = vcf::Record::builder()
        .set_chromosome("sq0".parse()?)
        .set_position(Position::try_from(1)?)
        .set_reference_bases("A".parse()?)
        .build()?;

    writer.write_vcf_record(&header, &string_map, &record)?;

    Ok(())
}
