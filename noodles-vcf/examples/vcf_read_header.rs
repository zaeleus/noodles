//! Prints the header of a VCF file.
//!
//! The header is defined to be the meta lines (`##` prefix) and header line (`#`` prefix).
//!
//! The result is similar or matches the output of `bcftools view --header-only --no-version
//! <src>`. bcftools may add a PASS FILTER to the meta if it is missing.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_vcf as vcf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");
    let mut reader = File::open(src).map(BufReader::new).map(vcf::Reader::new)?;

    let header = reader.read_header()?;
    print!("{}", header);

    Ok(())
}
