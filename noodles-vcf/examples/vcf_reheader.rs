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

use noodles_vcf::{
    self as vcf,
    header::{
        record::{Key, Value},
        Record,
    },
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");
    let mut reader = File::open(src).map(BufReader::new).map(vcf::Reader::new)?;

    let mut header: vcf::Header = reader.read_header()?.parse()?;

    let record = Record::new(
        Key::Other(String::from("comment")),
        Value::String(String::from("a comment added by noodles-vcf")),
    );

    header.insert(record);

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = vcf::Writer::new(handle);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;
        writer.write_record(&record)?;
    }

    Ok(())
}
