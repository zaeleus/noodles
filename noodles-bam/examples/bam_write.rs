//! Creates a new BAM file.
//!
//! This writes a SAM header and three unmapped records to stdout.
//!
//! Verify the output by piping to `samtools view -h --no-PG`.

use std::io;

use noodles_bam as bam;
use noodles_sam::{
    self as sam,
    header::{self, Program, ReferenceSequence},
};

fn main() -> io::Result<()> {
    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = bam::Writer::new(handle);

    let header = sam::Header::builder()
        .set_header(header::header::Header::default())
        .add_reference_sequence(ReferenceSequence::new(String::from("sq0"), 8))
        .add_reference_sequence(ReferenceSequence::new(String::from("sq1"), 13))
        .add_reference_sequence(ReferenceSequence::new(String::from("sq2"), 21))
        .add_program(Program::new(String::from("noodles-bam")))
        .add_comment("an example BAM written by noodles-bam")
        .build();

    writer.write_header(&header)?;

    for _ in 0..3 {
        let record = sam::Record::default();
        writer.write_record(header.reference_sequences(), &record)?;
    }

    Ok(())
}
