//! Counts the number of records in a VCF file.
//!
//! The result matches the output of `bcftools view --no-header | wc -l`.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_vcf as vcf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(|f| vcf::Reader::new(BufReader::new(f)))?;
    reader.read_header()?;

    let count = reader.records().count();
    println!("{}", count);

    Ok(())
}
