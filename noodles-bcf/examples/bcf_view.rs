//! Prints a BCF file in the VCF format.
//!
//! The result matches the output of `bcftools view <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufWriter},
};

use noodles_bcf as bcf;
use noodles_vcf as vcf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bcf::Reader::new)?;
    reader.read_file_format()?;
    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = vcf::Writer::new(BufWriter::new(stdout));

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    Ok(())
}
