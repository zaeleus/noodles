//! Reading records with only the spans of reference sequences that slices use.

use std::{error::Error, io};

use noodles_core::Position;
use noodles_cram as cram;
use noodles_fasta::{self as fasta, repository::adapters};
use noodles_sam::{
    self as sam,
    alignment::{
        RecordBuf,
        io::Write as _,
        record::{
            Flags,
            cigar::{
                Op,
                op::Kind::{self, Deletion, Insertion, Match, Skip, SoftClip},
            },
        },
    },
};

type References = [(&'static str, Vec<u8>)];

const INDELS: &[(Kind, usize)] = &[
    (SoftClip, 2),
    (Match, 10),
    (Insertion, 1),
    (Match, 5),
    (Deletion, 2),
    (Match, 8),
];

const SPLICED: &[(Kind, usize)] = &[(Match, 10), (Skip, 3), (Match, 10)];

#[test]
fn test_single_reference_slice() -> Result<(), Box<dyn Error>> {
    // `sq0` is long enough for the slice to read only its span.
    let mut references = [("sq0", reference_sequence(1000))];

    let records = [
        record("r0", &references, 0, 100, INDELS),
        record("r1", &references, 0, 110, SPLICED),
    ];

    let src = write(&references, &records)?;

    let indexed = indexed_repository(&references)?;
    assert_eq!(read(&src, indexed.clone())?, records);
    // `sq0` was not read whole, which would have cached it.
    assert!(indexed.is_empty());

    // A reference base that differs in the slice span fails the slice reference checksum.
    references[0].1[104] = b'N';

    for repository in [repository(&references), indexed_repository(&references)?] {
        let e = read(&src, repository).unwrap_err();
        assert_eq!(e.kind(), io::ErrorKind::InvalidData);
    }

    Ok(())
}

#[test]
fn test_multi_reference_slice() -> Result<(), Box<dyn Error>> {
    // `sq0` is long enough for the slice to read only the span of its records. `sq1` is not.
    let references = [
        ("sq0", reference_sequence(1000)),
        ("sq1", reference_sequence(50)),
    ];

    let records = [
        record("r0", &references, 0, 100, INDELS),
        record("r1", &references, 1, 10, &[(Match, 30)]),
        record("r2", &references, 1, 20, INDELS),
    ];

    let src = write(&references, &records)?;

    let indexed = indexed_repository(&references)?;
    assert_eq!(read(&src, indexed.clone())?, records);
    // Only `sq1` was read whole.
    assert_eq!(indexed.len(), 1);

    Ok(())
}

/// Returns a pseudorandom sequence, so that reading the wrong part of it changes the bases of
/// records.
fn reference_sequence(len: usize) -> Vec<u8> {
    let mut state = 1u32;

    (0..len)
        .map(|_| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            b"ACGT"[(state >> 16) as usize % 4]
        })
        .collect()
}

/// Returns a mapped record whose aligned bases are those of the reference sequence.
fn record(
    name: &str,
    references: &References,
    reference_sequence_id: usize,
    alignment_start: usize,
    ops: &[(Kind, usize)],
) -> RecordBuf {
    let reference_sequence = &references[reference_sequence_id].1;

    let mut sequence = Vec::new();
    let mut i = alignment_start - 1;

    for &(kind, len) in ops {
        match kind {
            Match => {
                sequence.extend(&reference_sequence[i..i + len]);
                i += len;
            }
            Deletion | Skip => i += len,
            _ => sequence.extend(b"N".repeat(len)),
        }
    }

    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::empty())
        .set_reference_sequence_id(reference_sequence_id)
        .set_alignment_start(Position::new(alignment_start).unwrap())
        .set_cigar(ops.iter().map(|&(kind, len)| Op::new(kind, len)).collect())
        .set_quality_scores(vec![30; sequence.len()].into())
        .set_sequence(sequence.into())
        .build()
}

fn fasta_records(references: &References) -> Vec<fasta::Record> {
    references
        .iter()
        .map(|(name, sequence)| {
            let definition = fasta::record::Definition::new(*name, None);
            fasta::Record::new(definition, sequence.clone().into())
        })
        .collect()
}

/// Returns a repository that reads whole sequences.
fn repository(references: &References) -> fasta::Repository {
    fasta::Repository::new(fasta_records(references))
}

/// Returns a repository of an indexed FASTA, which reads ranges of sequences.
fn indexed_repository(references: &References) -> Result<fasta::Repository, Box<dyn Error>> {
    let mut writer = fasta::io::Writer::new(Vec::new());

    for record in fasta_records(references) {
        writer.write_record(&record)?;
    }

    let src = writer.into_inner();

    let mut indexer = fasta::io::Indexer::new(&src[..]);
    let mut index = Vec::new();

    while let Some(record) = indexer.index_record()? {
        index.push(record);
    }

    let reader = fasta::io::IndexedReader::new(io::Cursor::new(src), index.into());
    Ok(fasta::Repository::new(adapters::IndexedReader::new(reader)))
}

fn write(references: &References, records: &[RecordBuf]) -> Result<Vec<u8>, Box<dyn Error>> {
    let header: sam::Header = references
        .iter()
        .map(|(name, sequence)| format!("@SQ\tSN:{name}\tLN:{}\n", sequence.len()))
        .collect::<String>()
        .parse()?;

    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository(references))
        .build_from_writer(Vec::new());

    writer.write_header(&header)?;

    for record in records {
        writer.write_alignment_record(&header, record)?;
    }

    writer.try_finish(&header)?;

    Ok(writer.into_inner())
}

fn read(src: &[u8], repository: fasta::Repository) -> io::Result<Vec<RecordBuf>> {
    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_reader(src);

    let header = reader.read_header()?;
    reader.records(&header).collect()
}
