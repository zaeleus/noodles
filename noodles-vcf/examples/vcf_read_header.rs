//! Prints the header of a VCF file.
//!
//! The header is defined to be the meta lines (`##` prefix) and header line (`#`` prefix).
//!
//! The result is similar to or matches the output of `bcftools head <src>`. bcftools may add a
//! PASS FILTER to the meta if it is missing.

use std::{env, io};

use noodles_vcf as vcf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = vcf::io::reader::Builder::default().build_from_path(src)?;

    let header = reader.read_header()?;
    print!("{header}");

    Ok(())
}
