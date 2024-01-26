//! Prints the header of a BCF file.
//!
//! The result matches the output of `bcftools head <src>`.

use std::{env, io};

use noodles_bcf as bcf;
use noodles_vcf as vcf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = bcf::io::reader::Builder::default().build_from_path(src)?;
    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = vcf::io::Writer::new(stdout);
    writer.write_header(&header)?;

    Ok(())
}
