//! Prints a BCF file in the VCF format.
//!
//! The result matches the output of `bcftools view --no-header <src>`.

use std::{env, fs::File};

use noodles_bcf as bcf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bcf::Reader::new)?;
    reader.read_file_format()?;

    let raw_header = reader.read_header()?;
    let header = raw_header.parse()?;
    let string_map = raw_header.parse()?;

    for result in reader.records() {
        let record = result?;
        let vcf_record = record.try_into_vcf_record(&header, &string_map)?;
        println!("{}", vcf_record);
    }

    Ok(())
}
