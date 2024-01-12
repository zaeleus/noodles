//! Prints the header of a BAM file.
//!
//! A BAM file header is a SAM header.
//!
//! The result is similar to the output of `samtools head <src>`.

use std::{env, io};

use noodles_bam as bam;
use noodles_sam as sam;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = bam::reader::Builder.build_from_path(src)?;
    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = sam::Writer::new(stdout);
    writer.write_header(&header)?;

    Ok(())
}
