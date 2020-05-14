//! Creates a new SAM file.
//!
//! This writes a SAM header and three unmapped records to stdout.
//!
//! Verify the output by piping to `samtools view -h --no-PG`.

use std::io;

use noodles_sam::{
    self as sam,
    header::{Program, ReferenceSequence},
};

fn main() -> io::Result<()> {
    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = sam::Writer::new(handle);

    let header = sam::Header::builder()
        .add_reference_sequence(ReferenceSequence::new(String::from("sq0"), 8))
        .add_reference_sequence(ReferenceSequence::new(String::from("sq1"), 13))
        .add_reference_sequence(ReferenceSequence::new(String::from("sq2"), 21))
        .add_program(Program::new(String::from("noodles-sam")))
        .add_comment("an example SAM written by noodles-sam")
        .build();

    writer.write_header(&header)?;

    for _ in 0..3 {
        let record = sam::Record::default();
        writer.write_record(&record)?;
    }

    Ok(())
}
