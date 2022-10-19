//! Replaces the SAM header of a BAM file.
//!
//! This is similar to the functionality of `samtools reheader`.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::{env, fs::File, io};

use noodles_bam as bam;
use noodles_sam as sam;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bam::Reader::new)?;
    let mut header: sam::Header = reader.read_header()?.parse()?;
    reader.read_reference_sequences()?;

    header.add_comment("a comment added by noodles-bam");

    let stdout = io::stdout().lock();
    let mut writer = bam::Writer::new(stdout);

    writer.write_header(&header)?;
    writer.write_reference_sequences(header.reference_sequences())?;

    for result in reader.records() {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    Ok(())
}
