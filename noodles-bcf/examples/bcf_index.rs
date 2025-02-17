//! Builds and writes a coordinate-sorted index (CSI) from a BCF file.
//!
//! This writes the output to stdout rather than `<src>.csi`.
//!
//! The output is similar to the output of `bcftools index <src>`.

use std::{
    env,
    io::{self, BufWriter},
};

use noodles_bcf as bcf;
use noodles_csi as csi;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let index = bcf::fs::index(src)?;

    let stdout = io::stdout().lock();
    let mut writer = csi::io::Writer::new(BufWriter::new(stdout));
    writer.write_index(&index)?;

    Ok(())
}
