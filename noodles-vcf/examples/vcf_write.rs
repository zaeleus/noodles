//! This writes a VCF header and one (sparse) record to stdout.
//!
//! Verify the output by piping to `bcftools view --no-version`.

use std::io;

use noodles_vcf::{
    self as vcf,
    header::record::value::{map::Contig, Map},
    record::Position,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let stdout = io::stdout().lock();
    let mut writer = vcf::io::Writer::new(stdout);

    let header = vcf::Header::builder()
        .add_contig("sq0", Map::<Contig>::new())
        .build();

    writer.write_header(&header)?;

    let record = vcf::Record::builder()
        .set_chromosome("sq0")
        .set_position(Position::from(1))
        .set_reference_bases("A".parse()?)
        .build()?;

    writer.write_record(&header, &record)?;

    Ok(())
}
