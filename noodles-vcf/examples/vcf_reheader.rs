//! Replaces the header of a VCF file.
//!
//! This is similar to the functionality of `bcftools reheader --header <header-src> <src>`.
//!
//! Verify the output by piping to `bcftools view --no-version --header`.

use std::{
    env,
    fs::File,
    io::{self, BufReader, BufWriter},
};

use noodles_vcf as vcf;

fn add_comment(header: &mut vcf::Header) -> Result<(), Box<dyn std::error::Error>> {
    use vcf::header::record::Value;

    header.insert(
        "comment".parse()?,
        Value::from("a comment added by noodles-vcf"),
    )?;

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(vcf::io::Reader::new)?;

    let mut header = reader.read_header()?;
    add_comment(&mut header)?;

    let stdout = io::stdout().lock();
    let mut writer = vcf::io::Writer::new(BufWriter::new(stdout));

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    Ok(())
}
