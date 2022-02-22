//! Creates a new CRAM file.
//!
//! This writes a file definition, a header container built from a SAM header, one mapped record,
//! and one unmapped record to stdout.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::io;

use md5::{Digest, Md5};

use noodles_bam as bam;
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

fn build_header(
    reference_sequence_records: &[fasta::Record],
) -> Result<sam::Header, Box<dyn std::error::Error>> {
    let mut builder = sam::Header::builder()
        .set_header(header::header::Header::default())
        .add_program(Program::new("noodles-cram"))
        .add_comment("an example CRAM written by noodles-cram");

    for record in reference_sequence_records {
        let sequence = record.sequence();

        let mut hasher = Md5::new();
        hasher.update(&sequence);
        let md5_checksum = Md5Checksum::from(<[u8; 16]>::from(hasher.finalize()));

        let reference_sequence = ReferenceSequence::builder()
            .set_name(record.name().parse()?)
            .set_length(sequence.len() as i32)
            .set_md5_checksum(md5_checksum)
            .build()?;

        builder = builder.add_reference_sequence(reference_sequence);
    }

    Ok(builder.build())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let reference_sequence_records: Vec<_> = fasta::Reader::new(FASTA_DATA)
        .records()
        .collect::<Result<_, _>>()?;

    let header = build_header(&reference_sequence_records)?;

    let stdout = io::stdout();
    let handle = stdout.lock();

    let mut writer = cram::Writer::new(handle, reference_sequence_records);
    writer.write_file_definition()?;
    writer.write_file_header(&header)?;

    let record = cram::Record::builder()
        .set_bam_flags(sam::record::Flags::empty())
        .set_flags(cram::record::Flags::QUALITY_SCORES_STORED_AS_ARRAY)
        .set_reference_sequence_id(bam::record::ReferenceSequenceId::from(0))
        .set_alignment_start(sam::record::Position::try_from(1)?)
        .set_read_length(4)
        .set_bases(b"TTCA".to_vec())
        .set_quality_scores(vec![45, 35, 43, 50])
        .build();

    writer.write_record(record)?;

    let record = cram::Record::default();
    writer.write_record(record)?;

    Ok(())
}
