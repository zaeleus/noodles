//! Creates a new BAM file.
//!
//! This writes a SAM header, reference sequences, and three unmapped records to stdout.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::io;

use noodles_bam as bam;
use noodles_sam::{
    self as sam,
    header::{self, Program, ReferenceSequence},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = bam::Writer::new(handle);

    let header = sam::Header::builder()
        .set_header(header::header::Header::default())
        .add_reference_sequence(ReferenceSequence::new("sq0".parse()?, 8)?)
        .add_reference_sequence(ReferenceSequence::new("sq1".parse()?, 13)?)
        .add_reference_sequence(ReferenceSequence::new("sq2".parse()?, 21)?)
        .add_program(Program::new("noodles-bam"))
        .add_comment("an example BAM written by noodles-bam")
        .build();

    writer.write_header(&header)?;
    writer.write_reference_sequences(header.reference_sequences())?;

    for _ in 0..3 {
        let record = bam::Record::default();
        writer.write_record(&record)?;
    }

    Ok(())
}
