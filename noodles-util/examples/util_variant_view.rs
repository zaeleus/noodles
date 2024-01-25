//! Prints a variant file in the VCF format.
//!
//! The result matches the output of `bcftools view <src>`.

use std::{
    env,
    io::{self, BufWriter},
};

use noodles_util::variant;
use noodles_vcf as vcf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let builder = variant::reader::Builder::default();

    let mut reader = if src == "-" {
        let stdin = io::stdin().lock();
        builder.build_from_reader(stdin)?
    } else {
        builder.build_from_path(src)?
    };

    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = vcf::io::Writer::new(BufWriter::new(stdout));

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    Ok(())
}
