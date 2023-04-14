//! Queries a BCF file with a given region.
//!
//! The input BCF must have an index in the same directory.
//!
//! The result matches the output of `bcftools view --no-header <src> <region>`.

use std::{
    env,
    fs::File,
    io::{self, BufWriter},
    path::PathBuf,
};

use noodles_bcf::{self as bcf, header::StringMaps};
use noodles_csi as csi;
use noodles_vcf as vcf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args();

    let src = args.nth(1).map(PathBuf::from).expect("missing src");
    let raw_region = args.next().expect("missing region");

    let mut reader = File::open(&src).map(bcf::Reader::new)?;
    reader.read_file_format()?;
    let raw_header = reader.read_header()?;

    let header: vcf::Header = raw_header.parse()?;
    let string_maps: StringMaps = raw_header.parse()?;

    let index = csi::read(src.with_extension("bcf.csi"))?;

    let region = raw_region.parse()?;
    let query = reader.query(string_maps.contigs(), &index, &region)?;

    let stdout = io::stdout().lock();
    let mut writer = vcf::Writer::new(BufWriter::new(stdout));

    for result in query {
        let record = result?;
        let vcf_record = record.try_into_vcf_record(&header, &string_maps)?;
        writer.write_record(&header, &vcf_record)?;
    }

    Ok(())
}
