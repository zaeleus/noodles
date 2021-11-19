//! Queries a VCF file with a given region.
//!
//! The input VCF must have an associated index in the same directory.
//!
//! The result matches the output `bcftools view --no-header <src> <region>`.

use std::{env, fs::File, path::PathBuf};

use noodles_bgzf as bgzf;
use noodles_tabix as tabix;
use noodles_vcf as vcf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args();

    let src = args.nth(1).map(PathBuf::from).expect("missing src");
    let raw_region = args.next().expect("missing region");

    let mut reader = File::open(&src)
        .map(bgzf::Reader::new)
        .map(vcf::Reader::new)?;

    let header = reader.read_header()?.parse()?;

    let index = tabix::read(src.with_extension("gz.tbi"))?;
    let region = raw_region.parse()?;

    let query = reader.query(&header, &index, &region)?;

    for result in query {
        let record = result?;
        println!("{}", record);
    }

    Ok(())
}
