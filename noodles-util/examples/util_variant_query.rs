//! Queries a variant file with the given region.
//!
//! The input must have an associated index in the same directory.
//!
//! The result matches the output of `bcftools view <src> <region>`.

use std::{
    env,
    io::{self, BufWriter},
};

use noodles_util::variant;
use noodles_vcf::{self as vcf, variant::io::Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let region = args.next().expect("missing region").parse()?;

    let mut reader = variant::io::indexed_reader::Builder::default().build_from_path(src)?;
    let header = reader.read_header()?;

    let query = reader.query(&header, &region)?;

    let stdout = io::stdout().lock();
    let mut writer = vcf::io::Writer::new(BufWriter::new(stdout));

    writer.write_header(&header)?;

    for result in query {
        let record = result?;
        writer.write_variant_record(&header, &record)?;
    }

    Ok(())
}
