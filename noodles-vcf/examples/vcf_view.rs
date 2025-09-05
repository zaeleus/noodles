//! Validates and prints a VCF file to stdout.
//!
//! The result matches the output of `bcftools view <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufReader, BufWriter},
};

use noodles_vcf as vcf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(vcf::io::Reader::new)?;

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
