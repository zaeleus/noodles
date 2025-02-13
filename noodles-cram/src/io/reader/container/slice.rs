mod header;

use std::io;

use bstr::BString;
use bytes::Bytes;
use noodles_core::Position;
use noodles_fasta as fasta;
use noodles_sam as sam;

use self::header::get_header;
use super::{read_block, Block};
use crate::{
    calculate_normalized_sequence_digest,
    container::{
        block::ContentType, compression_header::PreservationMap, slice::Header, CompressionHeader,
        ReferenceSequenceContext,
    },
    io::BitReader,
    record::{resolve, Feature},
    Record,
};

/// A container slice.
///
/// A slice contains a header, a core data block, and one or more external blocks. This is where
/// the CRAM records are stored.
pub struct Slice {
    header: Header,
    src: Bytes,
}

impl Slice {
    pub(crate) fn header(&self) -> &Header {
        &self.header
    }

    /// Reads and returns a list of raw records in this slice.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_cram::{self as cram, io::reader::Container};
    ///
    /// let data = [];
    /// let mut reader = cram::io::Reader::new(&data[..]);
    /// reader.read_header()?;
    ///
    /// let mut container = Container::default();
    ///
    /// while reader.read_container(&mut container)? != 0 {
    ///     let compression_header = container.compression_header()?;
    ///
    ///     for result in container.slices() {
    ///         let slice = result?;
    ///         let records = slice.records(&compression_header)?;
    ///         // ...
    ///     }
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records(&self, compression_header: &CompressionHeader) -> io::Result<Vec<Record>> {
        use crate::io::reader::record::ExternalDataReaders;

        let mut src = self.src.clone();

        let core_data_block = read_core_data_block(&mut src)?;

        let external_block_count = self.header.block_count() - 1;
        let external_blocks = read_external_blocks(&mut src, external_block_count)?;

        let core_data_reader = core_data_block
            .decode()
            .map(|buf| BitReader::new(Bytes::from(buf)))?;

        let mut external_data_readers = ExternalDataReaders::new();

        for block in &external_blocks {
            let reader = block.decode()?;
            external_data_readers.insert(block.content_id, Bytes::from(reader));
        }

        let mut record_reader = crate::io::reader::record::Reader::new(
            compression_header,
            core_data_reader,
            external_data_readers,
            self.header.reference_sequence_context(),
        );

        let record_count = self.header.record_count();

        let mut records = vec![Record::default(); record_count];

        let start_id = self.header.record_counter();
        let end_id = start_id + (record_count as u64);
        let ids = start_id..end_id;

        for (id, record) in ids.zip(&mut records) {
            record.id = id;
            record_reader.read_record(record)?;
        }

        Ok(records)
    }

    /// Resolves records.
    ///
    /// This resolves mates, read names, bases, and quality scores.
    pub fn resolve_records(
        &self,
        reference_sequence_repository: &fasta::Repository,
        header: &sam::Header,
        compression_header: &CompressionHeader,
        records: &mut [Record],
    ) -> io::Result<()> {
        let mut src = self.src.clone();

        read_core_data_block(&mut src)?;

        let external_block_count = self.header.block_count() - 1;
        let external_blocks = read_external_blocks(&mut src, external_block_count)?;

        resolve_mates(records)?;

        resolve_bases(
            reference_sequence_repository,
            header,
            compression_header,
            &self.header,
            &external_blocks,
            records,
        )?;

        resolve_quality_scores(records);

        Ok(())
    }
}

pub fn read_slice(src: &mut Bytes) -> io::Result<Slice> {
    let header = get_header(src)?;

    Ok(Slice {
        header,
        src: src.clone(),
    })
}

fn read_core_data_block(src: &mut Bytes) -> io::Result<Block> {
    let block = read_block(src)?;

    if block.content_type != ContentType::CoreData {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid block content type: expected {:?}, got {:?}",
                ContentType::CoreData,
                block.content_type
            ),
        ));
    }

    Ok(block)
}

fn read_external_blocks(src: &mut Bytes, len: usize) -> io::Result<Vec<Block>> {
    let mut external_blocks = Vec::with_capacity(len);

    for _ in 0..len {
        let block = read_block(src)?;

        if block.content_type != ContentType::ExternalData {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid block content type: expected {:?}, got {:?}",
                    ContentType::ExternalData,
                    block.content_type
                ),
            ));
        }

        external_blocks.push(block);
    }

    Ok(external_blocks)
}

fn resolve_mates(records: &mut [Record]) -> io::Result<()> {
    let mut mate_indices: Vec<_> = records
        .iter()
        .enumerate()
        .map(|(i, record)| record.distance_to_next_fragment().map(|len| i + len + 1))
        .collect();

    for i in 0..records.len() {
        let record = &mut records[i];

        if record.name().is_none() {
            // SAFETY: `u64::to_string` is always a valid read name.
            let name = BString::from(record.id().to_string());
            record.name = Some(name);
        }

        if mate_indices[i].is_none() {
            continue;
        }

        let mut j = i;

        while let Some(mate_index) = mate_indices[j] {
            let mid = j + 1;
            let (left, right) = records.split_at_mut(mid);

            let record = &mut left[j];
            let mate = &mut right[mate_index - mid];
            set_mate(record, mate);

            if mate.name().is_none() {
                mate.name = record.name().map(|name| name.into());
            }

            j = mate_index;
        }

        let (left, right) = records.split_at_mut(j);
        let record = &mut right[0];
        let mate = &mut left[i];
        set_mate(record, mate);

        // "The TLEN field is positive for the leftmost segment of the template, negative for the
        // rightmost, and the sign for any middle segment is undefined. If segments cover the same
        // coordinates then the choice of which is leftmost and rightmost is arbitrary..."
        let template_size = calculate_template_size(record, mate);
        records[i].template_length = template_size;

        let mut j = i;

        while let Some(mate_index) = mate_indices[j] {
            let record = &mut records[mate_index];
            record.template_length = -template_size;
            mate_indices[j] = None;
            j = mate_index;
        }
    }

    Ok(())
}

fn set_mate(record: &mut Record, mate: &mut Record) {
    set_mate_chunk(
        &mut record.bam_flags,
        &mut record.mate_reference_sequence_id,
        &mut record.mate_alignment_start,
        mate.bam_flags(),
        mate.reference_sequence_id(),
        mate.alignment_start(),
    );
}

fn calculate_template_size(record: &Record, mate: &Record) -> i32 {
    calculate_template_size_chunk(
        record.alignment_start(),
        record.read_length(),
        record.features(),
        mate.alignment_start(),
        mate.read_length(),
        mate.features(),
    )
}

fn set_mate_chunk(
    record_bam_bit_flags: &mut sam::alignment::record::Flags,
    record_next_fragment_reference_sequence_id: &mut Option<usize>,
    record_next_mate_alignment_start: &mut Option<Position>,
    mate_bam_bit_flags: sam::alignment::record::Flags,
    mate_reference_sequence_id: Option<usize>,
    mate_alignment_start: Option<Position>,
) {
    if mate_bam_bit_flags.is_reverse_complemented() {
        *record_bam_bit_flags |= sam::alignment::record::Flags::MATE_REVERSE_COMPLEMENTED;
    }

    if mate_bam_bit_flags.is_unmapped() {
        *record_bam_bit_flags |= sam::alignment::record::Flags::MATE_UNMAPPED;
    }

    *record_next_fragment_reference_sequence_id = mate_reference_sequence_id;
    *record_next_mate_alignment_start = mate_alignment_start;
}

// _Sequence Alignment/Map Format Specification_ (2021-06-03) § 1.4.9 "TLEN"
fn calculate_template_size_chunk(
    record_alignment_start: Option<Position>,
    record_read_length: usize,
    record_features: &[Feature],
    mate_alignment_start: Option<Position>,
    mate_read_length: usize,
    mate_features: &[Feature],
) -> i32 {
    use crate::record::calculate_alignment_span;

    fn alignment_end(
        alignment_start: Option<Position>,
        read_length: usize,
        features: &[Feature],
    ) -> Option<Position> {
        alignment_start.and_then(|start| {
            let span = calculate_alignment_span(read_length, features);
            let end = usize::from(start) + span - 1;
            Position::new(end)
        })
    }

    let Some(start) = record_alignment_start
        .min(mate_alignment_start)
        .map(usize::from)
    else {
        return 0;
    };

    let record_alignment_end =
        alignment_end(record_alignment_start, record_read_length, record_features);
    let mate_alignment_end = alignment_end(mate_alignment_start, mate_read_length, mate_features);

    let end = record_alignment_end
        .max(mate_alignment_end)
        .map(usize::from)
        .expect("invalid end position");

    // "...the absolute value of TLEN equals the distance between the mapped end of the template
    // and the mapped start of the template, inclusively..."
    if start > end {
        (start - end + 1) as i32
    } else {
        (end - start + 1) as i32
    }
}

fn resolve_bases(
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    compression_header: &CompressionHeader,
    slice_header: &Header,
    slice_external_data_blocks: &[Block],
    records: &mut [Record],
) -> io::Result<()> {
    let preservation_map = compression_header.preservation_map();

    let slice_reference_sequence = get_slice_reference_sequence(
        reference_sequence_repository,
        header,
        preservation_map,
        slice_header,
        slice_external_data_blocks,
    )?;

    let is_reference_required = preservation_map.is_reference_required();

    for record in records {
        if record.bam_flags().is_unmapped() || record.cram_flags().decode_sequence_as_unknown() {
            continue;
        }

        let mut alignment_start = record.alignment_start.expect("invalid alignment start");

        let reference_sequence = if is_reference_required {
            if let Some(SliceReferenceSequence::External(reference_sequence_id, sequence)) =
                &slice_reference_sequence
            {
                if record.reference_sequence_id() == Some(*reference_sequence_id) {
                    Some(sequence.clone())
                } else {
                    // An invalid state?
                    todo!();
                }
            } else {
                let reference_sequence_name = record
                    .reference_sequence(header.reference_sequences())
                    .transpose()?
                    .map(|(name, _)| name)
                    .expect("invalid reference sequence ID");

                let sequence = reference_sequence_repository
                    .get(reference_sequence_name)
                    .transpose()?
                    .expect("invalid reference sequence name");

                Some(sequence)
            }
        } else if let Some(SliceReferenceSequence::Embedded(reference_start, sequence)) =
            &slice_reference_sequence
        {
            let offset_start = usize::from(alignment_start) - usize::from(*reference_start) + 1;
            alignment_start = Position::try_from(offset_start)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Some(sequence.clone())
        } else {
            None
        };

        let substitution_matrix = preservation_map.substitution_matrix();

        resolve::resolve_bases(
            reference_sequence.as_ref(),
            substitution_matrix,
            &record.features,
            alignment_start,
            record.read_length(),
            &mut record.sequence,
        )?;
    }

    Ok(())
}

enum SliceReferenceSequence {
    External(usize, fasta::record::Sequence),
    Embedded(Position, fasta::record::Sequence),
}

fn get_slice_reference_sequence(
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    preservation_map: &PreservationMap,
    slice_header: &Header,
    slice_external_data_blocks: &[Block],
) -> io::Result<Option<SliceReferenceSequence>> {
    let reference_sequence_context = slice_header.reference_sequence_context();

    let ReferenceSequenceContext::Some(context) = reference_sequence_context else {
        return Ok(None);
    };

    let is_reference_required = preservation_map.is_reference_required();
    let embedded_reference_bases_block_content_id =
        slice_header.embedded_reference_bases_block_content_id();

    if is_reference_required {
        let reference_sequence_name = header
            .reference_sequences()
            .get_index(context.reference_sequence_id())
            .map(|(name, _)| name)
            .expect("invalid slice reference sequence ID");

        let sequence = reference_sequence_repository
            .get(reference_sequence_name)
            .transpose()?
            .expect("invalid slice reference sequence name");

        // § 11 "Reference sequences" (2021-11-15): "All CRAM reader implementations are
        // expected to check for reference MD5 checksums and report any missing or
        // mismatching entries."
        let start = context.alignment_start();
        let end = context.alignment_end();

        let actual_md5 = calculate_normalized_sequence_digest(&sequence[start..=end]);
        let expected_md5 = slice_header.reference_md5();

        if actual_md5 != expected_md5 {
            return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "reference sequence checksum mismatch: expected {expected_md5:?}, got {actual_md5:?}"
                        ),
                    ));
        }

        Ok(Some(SliceReferenceSequence::External(
            context.reference_sequence_id(),
            sequence,
        )))
    } else if let Some(block_content_id) = embedded_reference_bases_block_content_id {
        let block = slice_external_data_blocks
            .iter()
            .find(|block| block.content_id == block_content_id)
            .expect("invalid block content ID");

        let data = block.decode()?;
        let sequence = fasta::record::Sequence::from(data);

        let reference_start = context.alignment_start();
        Ok(Some(SliceReferenceSequence::Embedded(
            reference_start,
            sequence,
        )))
    } else {
        Ok(None)
    }
}

fn resolve_quality_scores(records: &mut [Record]) {
    for record in records {
        if !record.flags().is_unmapped()
            && !record.cram_flags().are_quality_scores_stored_as_array()
        {
            resolve::resolve_quality_scores(
                &record.features,
                record.read_length(),
                &mut record.quality_scores,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;

    use super::*;
    use crate::{
        calculate_normalized_sequence_digest, container::block::CompressionMethod, record::Flags,
    };

    #[test]
    fn test_resolve_mates() -> io::Result<()> {
        let mut records = vec![
            Record {
                id: 1,
                cram_flags: Flags::HAS_MATE_DOWNSTREAM,
                reference_sequence_id: Some(2),
                read_length: 4,
                alignment_start: Position::new(5),
                distance_to_mate: Some(0),
                ..Default::default()
            },
            Record {
                id: 2,
                cram_flags: Flags::HAS_MATE_DOWNSTREAM,
                reference_sequence_id: Some(2),
                read_length: 4,
                alignment_start: Position::new(8),
                distance_to_mate: Some(1),
                ..Default::default()
            },
            Record {
                id: 3,
                ..Default::default()
            },
            Record {
                id: 4,
                reference_sequence_id: Some(2),
                read_length: 4,
                alignment_start: Position::new(13),
                ..Default::default()
            },
        ];

        resolve_mates(&mut records)?;

        let name_1 = b"1".as_bstr();

        assert_eq!(records[0].name(), Some(name_1));
        assert_eq!(
            records[0].next_fragment_reference_sequence_id(),
            records[1].reference_sequence_id()
        );
        assert_eq!(
            records[0].mate_alignment_start(),
            records[1].alignment_start(),
        );
        assert_eq!(records[0].template_size(), 12);

        assert_eq!(records[1].name(), Some(name_1));
        assert_eq!(
            records[1].next_fragment_reference_sequence_id(),
            records[3].reference_sequence_id()
        );
        assert_eq!(
            records[1].mate_alignment_start(),
            records[3].alignment_start(),
        );
        assert_eq!(records[1].template_size(), -12);

        let name_3 = b"3".as_bstr();
        assert_eq!(records[2].name(), Some(name_3));

        assert_eq!(records[3].name(), Some(name_1));
        assert_eq!(
            records[3].next_fragment_reference_sequence_id(),
            records[0].reference_sequence_id()
        );
        assert_eq!(
            records[3].mate_alignment_start(),
            records[0].alignment_start(),
        );
        assert_eq!(records[3].template_size(), -12);

        Ok(())
    }

    #[test]
    fn test_calculate_template_size() {
        use sam::alignment::record::Flags;

        // --> -->
        let record = Record {
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            alignment_start: Position::new(200),
            read_length: 50,
            ..Default::default()
        };

        assert_eq!(calculate_template_size(&record, &mate), 150);
        assert_eq!(calculate_template_size(&mate, &record), 150);

        // --> <--
        // This is the example given in _Sequence Alignment/Map Format Specification_ (2021-06-03)
        // § 1.4.9 "TLEN" (footnote 14).
        let record = Record {
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            bam_flags: Flags::REVERSE_COMPLEMENTED,
            alignment_start: Position::new(200),
            read_length: 50,
            ..Default::default()
        };

        assert_eq!(calculate_template_size(&record, &mate), 150);
        assert_eq!(calculate_template_size(&mate, &record), 150);

        // <-- -->
        let record = Record {
            bam_flags: Flags::REVERSE_COMPLEMENTED,
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            alignment_start: Position::new(200),
            read_length: 50,
            ..Default::default()
        };

        assert_eq!(calculate_template_size(&record, &mate), 150);
        assert_eq!(calculate_template_size(&mate, &record), 150);

        // <-- <--
        let record = Record {
            bam_flags: Flags::REVERSE_COMPLEMENTED,
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            bam_flags: Flags::REVERSE_COMPLEMENTED,
            alignment_start: Position::new(200),
            read_length: 50,
            ..Default::default()
        };

        assert_eq!(calculate_template_size(&record, &mate), 150);
        assert_eq!(calculate_template_size(&mate, &record), 150);

        // No alignment start position.
        let record = Record::default();
        assert_eq!(calculate_template_size(&record, &record), 0);
    }

    #[test]
    fn test_resolve_bases() -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZeroUsize;

        use sam::{
            alignment::record_buf::Sequence,
            header::record::value::map::{self, Map},
        };

        use crate::container::block::ContentType;

        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        let start = Position::try_from(1)?;
        let end = Position::try_from(2)?;
        let sequence = fasta::record::Sequence::from(b"ACGT".to_vec());
        let reference_md5 = calculate_normalized_sequence_digest(&sequence[start..=end]);

        let reference_sequence_repository = fasta::Repository::new(vec![fasta::Record::new(
            fasta::record::Definition::new("sq0", None),
            sequence,
        )]);

        let header = sam::Header::builder()
            .add_reference_sequence("sq0", Map::<map::ReferenceSequence>::new(SQ0_LN))
            .build();

        let compression_header = CompressionHeader::default();

        let slice_header = Header::builder()
            .set_reference_sequence_context(ReferenceSequenceContext::some(0, start, end))
            .set_reference_md5(reference_md5)
            .build();

        let slice_external_data_blocks = vec![Block {
            compression_method: CompressionMethod::None,
            content_type: ContentType::ExternalData,
            content_id: 1,
            uncompressed_size: 0,
            src: Vec::new(),
        }];

        let mut records = [Record {
            id: 1,
            bam_flags: sam::alignment::record::Flags::default(),
            reference_sequence_id: Some(0),
            read_length: 2,
            alignment_start: Some(Position::MIN),
            features: vec![Feature::Bases {
                position: Position::MIN,
                bases: vec![b'A', b'C'],
            }],
            ..Default::default()
        }];

        resolve_bases(
            &reference_sequence_repository,
            &header,
            &compression_header,
            &slice_header,
            &slice_external_data_blocks,
            &mut records,
        )?;

        let actual: Vec<_> = records.into_iter().map(|r| r.sequence).collect();
        let expected = [Sequence::from(vec![b'A', b'C'])];
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_resolve_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        use sam::alignment::record_buf::QualityScores;

        let mut records = [
            Record {
                id: 1,
                bam_flags: sam::alignment::record::Flags::empty(),
                read_length: 2,
                features: vec![Feature::Scores {
                    position: Position::try_from(1)?,
                    quality_scores: vec![8, 13],
                }],
                ..Default::default()
            },
            Record {
                id: 2,
                ..Default::default()
            },
            Record {
                id: 3,
                cram_flags: Flags::QUALITY_SCORES_STORED_AS_ARRAY,
                read_length: 2,
                quality_scores: QualityScores::from(vec![21, 34]),
                ..Default::default()
            },
        ];

        resolve_quality_scores(&mut records);

        let actual: Vec<_> = records.into_iter().map(|r| r.quality_scores).collect();

        let expected = [
            QualityScores::from(vec![8, 13]),
            QualityScores::default(),
            QualityScores::from(vec![21, 34]),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
