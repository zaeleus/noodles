//! This writes a VCF header and one (sparse) record to stdout.
//!
//! Verify the output by piping to `bcftools view --no-version`.

use std::io;

use noodles_core::Position;
use noodles_vcf::{
    self as vcf,
    header::record::value::{Map, map::Contig},
    variant::io::Write,
};

fn main() -> io::Result<()> {
    let stdout = io::stdout().lock();
    let mut writer = vcf::io::Writer::new(stdout);

    let header = vcf::Header::builder()
        .add_contig("sq0", Map::<Contig>::new())
        .build();

    writer.write_header(&header)?;

    let record = vcf::variant::RecordBuf::builder()
        .set_reference_sequence_name("sq0")
        .set_variant_start(Position::MIN)
        .set_reference_bases("A")
        .build();

    writer.write_variant_record(&header, &record)?;

    Ok(())
}
