//! Prints the header of a VCF file.
//!
//! The header is defined to be the meta lines (`##` prefix) and header line (`#`` prefix).
//!
//! The result is similar to or matches the output of `bcftools head <src>`. bcftools may add a
//! PASS FILTER to the meta if it is missing.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_vcf as vcf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(vcf::io::Reader::new)?;

    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = vcf::io::Writer::new(stdout);
    writer.write_header(&header)?;

    Ok(())
}
