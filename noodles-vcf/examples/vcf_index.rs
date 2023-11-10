//! Builds and write a tabix index from a bgzipped VCF file.
//!
//! This writes the output to stdout rather than `<src>.tbi`.
//!
//! The output is similar to the output of `bcftools index --tbi <src>`.

use std::{
    env,
    io::{self, BufWriter},
};

use noodles_tabix as tabix;
use noodles_vcf as vcf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let index = vcf::index(src)?;

    let stdout = io::stdout().lock();
    let mut writer = tabix::Writer::new(BufWriter::new(stdout));
    writer.write_index(&index)?;

    Ok(())
}
