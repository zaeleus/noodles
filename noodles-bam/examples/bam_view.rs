//! Prints a BAM file in the SAM format.
//!
//! The result matches the output of `samtools view <src>`.

use std::{env, io};

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::io::Write};

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = bam::reader::Builder.build_from_path(src)?;
    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = sam::Writer::new(stdout);

    for result in reader.records(&header) {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    writer.finish(&header)?;

    Ok(())
}
