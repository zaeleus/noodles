//! Prints a variant file in the VCF format.
//!
//! The result matches the output of `bcftools view <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufReader, BufWriter},
};

use noodles_util::variant;
use noodles_vcf::{self as vcf, variant::io::Write};

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = if src == "-" {
        variant::io::Reader::new(io::stdin().lock())?
    } else {
        File::open(src)
            .map(BufReader::new)
            .and_then(variant::io::Reader::new)?
    };

    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = vcf::io::Writer::new(BufWriter::new(stdout));

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;
        writer.write_variant_record(&header, record.as_ref())?;
    }

    Ok(())
}
