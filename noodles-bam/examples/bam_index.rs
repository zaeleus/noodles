//! Builds and writes a BAM index from a BAM file.
//!
//! The input BAM must be coordinate-sorted, i.e., `SO:coordinate`.
//!
//! This writes the output to stdout rather than `<src>.bai`.
//!
//! The output is similar to the output of `samtools index <src>`.

use std::{env, io};

use noodles_bam::{self as bam, bai};

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let index = bam::fs::index(src)?;

    let stdout = io::stdout().lock();
    let mut writer = bai::io::Writer::new(stdout);
    writer.write_index(&index)?;

    Ok(())
}
