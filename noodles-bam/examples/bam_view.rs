//! Prints a BAM file in the SAM format.
//!
//! The result matches the output of `samtools view <src>`.

use std::{env, fs::File, io};

use noodles_bam as bam;
use noodles_sam::{self as sam, AlignmentWriter};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bam::Reader::new)?;
    let header: sam::Header = reader.read_header()?.parse()?;
    reader.read_reference_sequences()?;

    let stdout = io::stdout().lock();
    let mut writer = sam::Writer::new(stdout);

    for result in reader.records() {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    writer.finish(&header)?;

    Ok(())
}
