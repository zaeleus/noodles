//! Creates a new CRAM file.
//!
//! This writes a file definition, a header container built from a SAM header, one mapped record,
//! and one unmapped record to stdout.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::{io, num::NonZeroUsize};

use noodles_core::Position;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::{
    self as sam,
    alignment::{
        io::Write,
        record_buf::{QualityScores, Sequence},
        RecordBuf,
    },
    header::record::value::{
        map::{Program, ReferenceSequence},
        Map,
    },
    record::cigar::{op::Kind, Op},
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
        .set_header(Default::default())
        .add_program("noodles-cram", Map::<Program>::default())
        .add_comment("an example CRAM written by noodles-cram");

    for record in reference_sequence_records {
        let name = record.name().parse()?;
        let length = NonZeroUsize::try_from(record.sequence().len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        let reference_sequence = Map::<ReferenceSequence>::new(length);
        builder = builder.add_reference_sequence(name, reference_sequence);
    }

    Ok(builder.build())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let reference_sequences = build_reference_sequences();

    let stdout = io::stdout().lock();

    let header = build_header(&reference_sequences)?;

    let repository = fasta::Repository::new(reference_sequences);
    let mut writer = cram::writer::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_with_writer(stdout);

    writer.write_header(&header)?;

    let record = RecordBuf::builder()
        .set_flags(sam::record::Flags::empty())
        .set_reference_sequence_id(1)
        .set_alignment_start(Position::MIN)
        .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
        .set_sequence(Sequence::from(b"TTCA".to_vec()))
        .set_quality_scores(QualityScores::from(vec![45, 35, 43, 50]))
        .build();

    writer.write_alignment_record(&header, &record)?;

    let record = RecordBuf::default();
    writer.write_alignment_record(&header, &record)?;

    writer.try_finish(&header)?;

    Ok(())
}
