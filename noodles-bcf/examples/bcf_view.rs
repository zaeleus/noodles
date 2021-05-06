//! Prints a BCF file in the VCF format.
//!
//! The result matches the output of `bcftools view --no-header <src>`.

use std::{env, fs::File};

use noodles_bcf as bcf;
use noodles_vcf as vcf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bcf::Reader::new)?;
    reader.read_file_format()?;

    let raw_header = reader.read_header()?;
    let header: vcf::Header = raw_header.parse()?;
    let string_map: bcf::header::StringMap = raw_header.parse()?;

    for result in reader.records() {
        let record = result.and_then(|r| r.try_into_vcf_record(&header, &string_map))?;
        println!("{}", record);
    }

    Ok(())
}
