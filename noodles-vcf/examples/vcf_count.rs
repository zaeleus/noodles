//! Counts the number of records in a VCF file.
//!
//! The result matches the output of `bcftools view --no-header <src> | wc -l`.

use std::env;

use noodles_vcf as vcf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = vcf::io::reader::Builder::default().build_from_path(src)?;
    reader.read_header()?;

    let mut n = 0;

    for result in reader.records() {
        let _ = result?;
        n += 1;
    }

    println!("{n}");

    Ok(())
}
