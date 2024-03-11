//! Replaces the header of a VCF file.
//!
//! This is similar to the functionality of `bcftools reheader --header <header-src> <src>`.
//!
//! Verify the output by piping to `bcftools view --no-version --header`.

use std::{
    env,
    io::{self, BufWriter},
};

use noodles_vcf::{self as vcf, variant::io::Write};

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

    let mut reader = vcf::io::reader::Builder::default().build_from_path(src)?;

    let mut header = reader.read_header()?;
    add_comment(&mut header)?;

    let stdout = io::stdout().lock();
    let mut writer = vcf::io::Writer::new(BufWriter::new(stdout));

    writer.write_header(&header)?;

    for result in reader.record_bufs(&header) {
        let record = result?;
        writer.write_variant_record(&header, &record)?;
    }

    Ok(())
}
