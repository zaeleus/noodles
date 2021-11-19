//! Validates and prints a VCF file to stdout.
//!
//! The result matches the output of `bcftools view <src>`.

use std::{env, fs::File, io::BufReader};

use noodles_vcf as vcf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(BufReader::new).map(vcf::Reader::new)?;
    let header = reader.read_header()?.parse()?;

    print!("{}", header);

    for result in reader.records(&header) {
        let record = result?;
        println!("{}", record);
    }

    Ok(())
}
