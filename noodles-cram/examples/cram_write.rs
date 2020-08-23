//! Creates a new CRAM file.
//!
//! This writes a file definition, a header container built from a SAM header, and one unmapped
//! record to stdout.
//!
//! Verify the output by piping to `samtools view -h --no-PG`.

use std::io;

use md5::{Digest, Md5};

use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::{
    self as sam,
    header::{self, reference_sequence::Md5Checksum, Program, ReferenceSequence},
};

static FASTA_DATA: &[u8] = b"\
>sq0
TTCACCCA
>sq1
GATCTTACTTTTT
>sq2
GGCGCCCCGCTGTGCAAAAAT
";

fn build_header(reference_sequence_records: &[fasta::Record]) -> sam::Header {
    let mut builder = sam::Header::builder()
        .set_header(header::header::Header::default())
        .add_program(Program::new(String::from("noodles-cram")))
        .add_comment("an example CRAM written by noodles-cram");

    for record in reference_sequence_records {
        let sequence = record.sequence();

        let mut hasher = Md5::new();
        hasher.update(&sequence);
        let md5_checksum = Md5Checksum::from(<[u8; 16]>::from(hasher.finalize()));

        let reference_sequence = ReferenceSequence::builder()
            .set_name(record.reference_sequence_name())
            .set_length(sequence.len() as i32)
            .set_md5_checksum(md5_checksum)
            .build();

        builder = builder.add_reference_sequence(reference_sequence);
    }

    builder.build()
}

fn main() -> io::Result<()> {
    let stdout = io::stdout();
    let handle = stdout.lock();

    let reference_sequence_records: Vec<_> = fasta::Reader::new(FASTA_DATA)
        .records()
        .collect::<Result<_, _>>()?;

    let mut writer = cram::Writer::new(handle, Vec::new());
    writer.write_file_definition()?;

    let header = build_header(&reference_sequence_records);
    writer.write_file_header(&header)?;

    let record = cram::Record::default();
    writer.write_record(&Vec::new(), record)?;

    Ok(())
}
