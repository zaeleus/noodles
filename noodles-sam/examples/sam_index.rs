//! Builds and writes a SAM index from a SAM file.
//!
//! The input SAM must be bgzip-compressed and coordinate-sorted, i.e., `SO:coordinate`.
//!
//! This writes the output to stdout rather than `<src>.csi`.
//!
//! The output is similar to the output of `samtools index -c <src>`.

use std::{env, io};

use noodles_csi as csi;
use noodles_sam as sam;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let index = sam::fs::index(src)?;

    let stdout = io::stdout().lock();
    let mut writer = csi::io::Writer::new(stdout);

    writer.write_index(&index)?;

    Ok(())
}
