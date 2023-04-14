//! Creates a new BCF file.
//!
//! This writes a BCF file format, VCF header, and a single VCF to stdout.
//!
//! Verify the output by piping to `bcftools view --no-version`.

use std::io;

use noodles_bcf::{self as bcf, header::StringMaps};
use noodles_vcf::{
    self as vcf,
    header::record::value::{
        map::{Contig, Filter},
        Map,
    },
    record::Position,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let stdout = io::stdout().lock();
    let mut writer = bcf::Writer::new(stdout);
    writer.write_file_format()?;

    let header = vcf::Header::builder()
        .add_filter("PASS", Map::<Filter>::pass())
        .add_contig("sq0".parse()?, Map::<Contig>::new())
        .build();

    writer.write_header(&header)?;

    let string_maps = StringMaps::from(&header);

    let record = vcf::Record::builder()
        .set_chromosome("sq0".parse()?)
        .set_position(Position::from(1))
        .set_reference_bases("A".parse()?)
        .build()?;

    writer.write_record(&header, &string_maps, &record)?;

    Ok(())
}
