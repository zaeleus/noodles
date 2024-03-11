//! Queries a VCF file with a given region.
//!
//! The input VCF must have an associated index in the same directory.
//!
//! The result matches the output `bcftools view --no-header <src> <region>`.

use std::{
    env,
    io::{self, BufWriter},
};

use noodles_vcf::{self as vcf, variant::io::Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args();

    let src = args.nth(1).expect("missing src");
    let raw_region = args.next().expect("missing region");

    let mut reader = vcf::io::indexed_reader::Builder::default().build_from_path(src)?;
    let header = reader.read_header()?;

    let region = raw_region.parse()?;
    let query = reader.query(&header, &region)?;

    let stdout = io::stdout().lock();
    let mut writer = vcf::io::Writer::new(BufWriter::new(stdout));

    for result in query {
        let record = result?;
        writer.write_variant_record(&header, &record)?;
    }

    Ok(())
}
