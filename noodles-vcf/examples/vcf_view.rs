//! Validates and prints a VCF file to stdout.
//!
//! The result matches the output of `bcftools view <src>`.

use std::{
    env,
    io::{self, BufWriter},
};

use noodles_vcf as vcf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = vcf::io::reader::Builder::default().build_from_path(src)?;
    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = vcf::io::Writer::new(BufWriter::new(stdout));

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    Ok(())
}
