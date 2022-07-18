//! Replaces the header of a VCF file.
//!
//! This is similar to the functionality of `bcftools reheader --header <header-src> <src>`.
//!
//! Verify the output by piping to `bcftools view --no-version --header`.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_vcf as vcf;

fn add_comment(header: &mut vcf::Header) {
    use vcf::header::record::{value, Key};

    header.insert(
        Key::from("comment"),
        value::Other::from("a comment added by noodles-vcf"),
    );
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(BufReader::new).map(vcf::Reader::new)?;

    let mut header = reader.read_header()?.parse()?;
    add_comment(&mut header);

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = vcf::Writer::new(handle);

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;
        writer.write_record(&record)?;
    }

    Ok(())
}
