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

fn build_header() -> sam::Header {
    use noodles_sam::header::reference_sequence::Tag;

    let builder = sam::Header::builder()
        .set_header(header::header::Header::default())
        .add_program(Program::new(String::from("noodles-cram")))
        .add_comment("an example CRAM written by noodles-cram");

    // TTCACCCA
    let mut sq0 = ReferenceSequence::new(String::from("sq0"), 8);
    sq0.insert(
        Tag::Md5Checksum,
        String::from("be19336b7e15968f7ac7dc82493d9cd8"),
    );
    let builder = builder.add_reference_sequence(sq0);

    // GATCTTACTTTTT
    let mut sq1 = ReferenceSequence::new(String::from("sq1"), 13);
    sq1.insert(
        Tag::Md5Checksum,
        String::from("d80f22a19aeeb623b3e4f746c762f21d"),
    );
    let builder = builder.add_reference_sequence(sq1);

    // GGCGCCCCGCTGTGCAAAAAT
    let mut sq2 = ReferenceSequence::new(String::from("sq2"), 21);
    sq2.insert(
        Tag::Md5Checksum,
        String::from("b00c61dfed4a92fdfb244d35790556eb"),
    );
    let builder = builder.add_reference_sequence(sq2);

    builder.build()
}

fn main() -> io::Result<()> {
    let stdout = io::stdout();
    let handle = stdout.lock();

    let mut writer = cram::Writer::new(handle);
    writer.write_file_definition()?;

    let header = build_header();
    writer.write_file_header(&header)?;

    Ok(())
}
