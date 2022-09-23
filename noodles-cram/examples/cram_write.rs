//! Creates a new CRAM file.
//!
//! This writes a file definition, a header container built from a SAM header, one mapped record,
//! and one unmapped record to stdout.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::io;

use md5::{Digest, Md5};

use noodles_core::Position;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::{
    self as sam,
    header::record::value::{
        map::{reference_sequence::Md5Checksum, Program, ReferenceSequence},
        Map,
    },
};
use sam::{alignment::Record, AlignmentWriter};

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
        .set_header(Default::default())
        .add_program(Map::<Program>::new("noodles-cram"))
        .add_comment("an example CRAM written by noodles-cram");

    for record in reference_sequence_records {
        let sequence = record.sequence();

        let name = record.name().parse()?;
        let len = sequence.len();

        let mut hasher = Md5::new();
        hasher.update(&sequence);
        let md5_checksum = Md5Checksum::from(<[u8; 16]>::from(hasher.finalize()));

        let reference_sequence = Map::<ReferenceSequence>::builder()
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
    let mut writer = cram::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_with_writer(handle);

    writer.write_file_definition()?;
    writer.write_file_header(&header)?;

    let record = Record::builder()
        .set_flags(sam::record::Flags::empty())
        .set_reference_sequence_id(1)
        .set_alignment_start(Position::MIN)
        .set_cigar("4M".parse()?)
        .set_sequence("TTCA".parse()?)
        .set_quality_scores("NDLS".parse()?)
        .build();

    writer.write_alignment_record(&header, &record)?;

    let record = Record::default();
    writer.write_alignment_record(&header, &record)?;

    writer.try_finish(&header)?;

    Ok(())
}
