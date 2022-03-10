//! Creates a new CRAM file.
//!
//! This writes a file definition, a header container built from a SAM header, one mapped record,
//! and one unmapped record to stdout.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::io;

use md5::{Digest, Md5};

use noodles_bam as bam;
use noodles_core::Position;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::{
    self as sam,
    header::{self, reference_sequence::Md5Checksum, Program, ReferenceSequence},
};

fn build_reference_sequences() -> Vec<fasta::Record> {
    use fasta::record::{Definition, Sequence};

    vec![
        fasta::Record::new(
            Definition::new("sq0", None),
            Sequence::from(b"TTCACCCA".to_vec()),
        ),
        fasta::Record::new(
            Definition::new("sq1", None),
            Sequence::from(b"GATCTTACTTTTT".to_vec()),
        ),
        fasta::Record::new(
            Definition::new("sq2", None),
            Sequence::from(b"GGCGCCCCGCTGTGCAAAAAT".to_vec()),
        ),
    ]
}

fn build_header(
    reference_sequence_records: &[fasta::Record],
) -> Result<sam::Header, Box<dyn std::error::Error>> {
    let mut builder = sam::Header::builder()
        .set_header(header::header::Header::default())
        .add_program(Program::new("noodles-cram"))
        .add_comment("an example CRAM written by noodles-cram");

    for record in reference_sequence_records {
        let sequence = record.sequence();

        let name = record.name().parse()?;
        let len = i32::try_from(sequence.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let mut hasher = Md5::new();
        hasher.update(&sequence);
        let md5_checksum = Md5Checksum::from(<[u8; 16]>::from(hasher.finalize()));

        let reference_sequence = ReferenceSequence::builder()
            .set_name(name)
            .set_length(len)
            .set_md5_checksum(md5_checksum)
            .build()?;

        builder = builder.add_reference_sequence(reference_sequence);
    }

    Ok(builder.build())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let reference_sequences = build_reference_sequences();

    let stdout = io::stdout();
    let handle = stdout.lock();

    let header = build_header(&reference_sequences)?;
    let repository = fasta::Repository::new(reference_sequences);
    let mut writer = cram::Writer::new(handle, repository, &header);

    writer.write_file_definition()?;
    writer.write_file_header(&header)?;

    let record = cram::Record::builder()
        .set_bam_flags(sam::record::Flags::empty())
        .set_flags(cram::record::Flags::QUALITY_SCORES_STORED_AS_ARRAY)
        .set_reference_sequence_id(bam::record::ReferenceSequenceId::from(0))
        .set_alignment_start(Position::try_from(1)?)
        .set_read_length(4)
        .set_bases("TTCA".parse()?)
        .set_quality_scores("NDLS".parse()?)
        .build();

    writer.write_record(record)?;

    let record = cram::Record::default();
    writer.write_record(record)?;

    Ok(())
}
