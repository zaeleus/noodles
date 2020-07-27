//! Creates a new CRAM file.
//!
//! This writes a file definition and a header container built from a SAM header to stdout.
//!
//! Verify the output by piping to `samtools view -h --no-PG`.

use std::io;

use noodles_cram as cram;
use noodles_sam::{
    self as sam,
    header::{self, Program, ReferenceSequence},
};

fn main() -> io::Result<()> {
    let stdout = io::stdout();
    let handle = stdout.lock();

    let mut writer = cram::Writer::new(handle);
    writer.write_file_definition()?;

    let header = sam::Header::builder()
        .set_header(header::header::Header::default())
        .add_reference_sequence(ReferenceSequence::new(String::from("sq0"), 8))
        .add_reference_sequence(ReferenceSequence::new(String::from("sq1"), 13))
        .add_reference_sequence(ReferenceSequence::new(String::from("sq2"), 21))
        .add_program(Program::new(String::from("noodles-cram")))
        .add_comment("an example CRAM written by noodles-cram")
        .build();

    writer.write_file_header(&header)?;

    Ok(())
}
