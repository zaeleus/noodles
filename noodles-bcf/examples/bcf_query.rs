//! Queries a BCF file with a given region.
//!
//! The input BCF must have an index in the same directory.
//!
//! The result matches the output of `bcftools view --no-header <src> <region>`.

use std::{
    env,
    fs::File,
    io::{self, BufWriter},
};

use noodles_bcf as bcf;
use noodles_vcf::{self as vcf, variant::io::Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let raw_region = args.next().expect("missing region");

    let mut reader = File::open(&src).map(bcf::io::Reader::new)?;
    let header = reader.read_header()?;

    let index = bcf::fs::read_associated_index(&src)?;
    let region = raw_region.parse()?;
    let query = reader.query(&header, &index, &region)?;

    let stdout = io::stdout().lock();
    let mut writer = vcf::io::Writer::new(BufWriter::new(stdout));

    for result in query.records() {
        let record = result?;
        writer.write_variant_record(&header, &record)?;
    }

    Ok(())
}
