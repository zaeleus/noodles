pub(crate) mod builder;
pub(crate) mod header;

pub use self::{builder::Builder, header::Header};

use std::io;

use noodles_core::Position;
use noodles_fasta as fasta;
use noodles_sam::{self as sam, alignment::record_buf::Name};

use super::{CompressionHeader, ReferenceSequenceContext};
use crate::{
    container::Block,
    io::BitReader,
    record::{resolve, Features},
    Record,
};

/// A CRAM data container slice.
///
/// A slice contains a header, a core data block, and one or more external blocks. This is where
/// the CRAM records are stored.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Slice {
    header: Header,
    core_data_block: Block,
    external_blocks: Vec<Block>,
}

impl Slice {
    pub(crate) fn builder() -> Builder {
        Builder::default()
    }

    pub(crate) fn new(header: Header, core_data_block: Block, external_blocks: Vec<Block>) -> Self {
        Self {
            header,
            core_data_block,
            external_blocks,
        }
    }

    pub(crate) fn header(&self) -> &Header {
        &self.header
    }

    pub(crate) fn core_data_block(&self) -> &Block {
        &self.core_data_block
    }

    pub(crate) fn external_blocks(&self) -> &[Block] {
        &self.external_blocks
    }

    /// Reads and returns a list of raw records in this slice.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_cram as cram;
    ///
    /// let data = [];
    /// let mut reader = cram::io::Reader::new(&data[..]);
    /// reader.read_header()?;
    ///
    /// while let Some(container) = reader.read_data_container()? {
    ///     for slice in container.slices() {
    ///         let records = slice.records(container.compression_header())?;
    ///         // ...
    ///     }
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records(&self, compression_header: &CompressionHeader) -> io::Result<Vec<Record>> {
        use crate::io::reader::record::ExternalDataReaders;

        let core_data_reader = self
            .core_data_block
            .decompressed_data()
            .map(BitReader::new)?;

        let mut external_data_readers = ExternalDataReaders::new();

        for block in self.external_blocks() {
            let reader = block.decompressed_data()?;
            external_data_readers.insert(block.content_id(), reader);
        }

        let mut record_reader = crate::io::reader::record::Reader::new(
            compression_header,
            core_data_reader,
            external_data_readers,
            self.header.reference_sequence_context(),
        );

        let record_count = self.header().record_count();

        let mut records = vec![Record::default(); record_count];

        let start_id = self.header().record_counter();
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
        resolve_mates(records)?;

        resolve_bases(
            reference_sequence_repository,
            header,
            compression_header,
            self,
            records,
        )?;

        resolve_quality_scores(records);

        Ok(())
    }
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
            let name = Name::from(record.id().to_string().into_bytes());
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
                mate.name = record.name().cloned();
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
        records[i].template_size = template_size;

        let mut j = i;

        while let Some(mate_index) = mate_indices[j] {
            let record = &mut records[mate_index];
            record.template_size = -template_size;
            mate_indices[j] = None;
            j = mate_index;
        }
    }

    Ok(())
}

fn set_mate(record: &mut Record, mate: &mut Record) {
    set_mate_chunk(
        &mut record.bam_bit_flags,
        &mut record.next_fragment_reference_sequence_id,
        &mut record.next_mate_alignment_start,
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
    record_bam_bit_flags: &mut sam::alignment::record_buf::Flags,
    record_next_fragment_reference_sequence_id: &mut Option<usize>,
    record_next_mate_alignment_start: &mut Option<Position>,
    mate_bam_bit_flags: sam::alignment::record_buf::Flags,
    mate_reference_sequence_id: Option<usize>,
    mate_alignment_start: Option<Position>,
) {
    if mate_bam_bit_flags.is_reverse_complemented() {
        *record_bam_bit_flags |= sam::alignment::record_buf::Flags::MATE_REVERSE_COMPLEMENTED;
    }

    if mate_bam_bit_flags.is_unmapped() {
        *record_bam_bit_flags |= sam::alignment::record_buf::Flags::MATE_UNMAPPED;
    }

    *record_next_fragment_reference_sequence_id = mate_reference_sequence_id;
    *record_next_mate_alignment_start = mate_alignment_start;
}

// _Sequence Alignment/Map Format Specification_ (2021-06-03) § 1.4.9 "TLEN"
fn calculate_template_size_chunk(
    record_alignment_start: Option<Position>,
    record_read_length: usize,
    record_features: &Features,
    mate_alignment_start: Option<Position>,
    mate_read_length: usize,
    mate_features: &Features,
) -> i32 {
    use crate::record::calculate_alignment_span;

    fn alignment_end(
        alignment_start: Option<Position>,
        read_length: usize,
        features: &Features,
    ) -> Option<Position> {
        alignment_start.and_then(|start| {
            let span = calculate_alignment_span(read_length, features);
            let end = usize::from(start) + span - 1;
            Position::new(end)
        })
    }

    let start = record_alignment_start
        .min(mate_alignment_start)
        .map(usize::from)
        .expect("invalid start positions");

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
    slice: &Slice,
    records: &mut [Record],
) -> io::Result<()> {
    enum SliceReferenceSequence {
        External(usize, fasta::record::Sequence),
        Embedded(usize, fasta::record::Sequence),
    }

    let preservation_map = compression_header.preservation_map();
    let is_reference_required = preservation_map.is_reference_required();
    let substitution_matrix = preservation_map.substitution_matrix();

    let slice_reference_sequence = if let ReferenceSequenceContext::Some(context) =
        slice.header().reference_sequence_context()
    {
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

            let actual_md5 = builder::calculate_normalized_sequence_digest(&sequence[start..=end]);
            let expected_md5 = slice.header().reference_md5();

            if actual_md5 != expected_md5 {
                return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "reference sequence checksum mismatch: expected {expected_md5:?}, got {actual_md5:?}"
                        ),
                    ));
            }

            Some(SliceReferenceSequence::External(
                context.reference_sequence_id(),
                sequence,
            ))
        } else if let Some(block_content_id) =
            slice.header().embedded_reference_bases_block_content_id()
        {
            let block = slice
                .external_blocks()
                .iter()
                .find(|block| block.content_id() == block_content_id)
                .expect("invalid block content ID");

            let data = block.decompressed_data()?;
            let sequence = fasta::record::Sequence::from(data);

            let offset = usize::from(context.alignment_start());
            Some(SliceReferenceSequence::Embedded(offset, sequence))
        } else {
            None
        }
    } else {
        None
    };

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
        } else if let Some(SliceReferenceSequence::Embedded(offset, sequence)) =
            &slice_reference_sequence
        {
            let start = usize::from(alignment_start) - offset + 1;
            alignment_start = Position::try_from(start)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Some(sequence.clone())
        } else {
            None
        };

        resolve::resolve_bases(
            reference_sequence.as_ref(),
            substitution_matrix,
            &record.features,
            alignment_start,
            record.read_length(),
            &mut record.bases,
        )?;
    }

    Ok(())
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
    use super::*;
    use crate::record::Flags;

    #[test]
    fn test_resolve_mates() -> Result<(), Box<dyn std::error::Error>> {
        let mut records = vec![
            Record::builder()
                .set_id(1)
                .set_flags(Flags::HAS_MATE_DOWNSTREAM)
                .set_reference_sequence_id(2)
                .set_read_length(4)
                .set_alignment_start(Position::try_from(5)?)
                .set_distance_to_next_fragment(0)
                .build(),
            Record::builder()
                .set_id(2)
                .set_flags(Flags::HAS_MATE_DOWNSTREAM)
                .set_reference_sequence_id(2)
                .set_read_length(4)
                .set_alignment_start(Position::try_from(8)?)
                .set_distance_to_next_fragment(1)
                .build(),
            Record::builder().set_id(3).build(),
            Record::builder()
                .set_id(4)
                .set_reference_sequence_id(2)
                .set_read_length(4)
                .set_alignment_start(Position::try_from(13)?)
                .build(),
        ];

        resolve_mates(&mut records)?;

        let name_1 = Name::from(b"1");

        assert_eq!(records[0].name(), Some(&name_1));
        assert_eq!(
            records[0].next_fragment_reference_sequence_id(),
            records[1].reference_sequence_id()
        );
        assert_eq!(
            records[0].mate_alignment_start(),
            records[1].alignment_start(),
        );
        assert_eq!(records[0].template_size(), 12);

        assert_eq!(records[1].name(), Some(&name_1));
        assert_eq!(
            records[1].next_fragment_reference_sequence_id(),
            records[3].reference_sequence_id()
        );
        assert_eq!(
            records[1].mate_alignment_start(),
            records[3].alignment_start(),
        );
        assert_eq!(records[1].template_size(), -12);

        let name_3 = Name::from(b"3");
        assert_eq!(records[2].name(), Some(&name_3));

        assert_eq!(records[3].name(), Some(&name_1));
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
    fn test_calculate_template_size() -> Result<(), noodles_core::position::TryFromIntError> {
        use sam::alignment::record_buf::Flags;

        // --> -->
        let record = Record::builder()
            .set_alignment_start(Position::try_from(100)?)
            .set_read_length(50)
            .build();

        let mate = Record::builder()
            .set_alignment_start(Position::try_from(200)?)
            .set_read_length(50)
            .build();

        assert_eq!(calculate_template_size(&record, &mate), 150);
        assert_eq!(calculate_template_size(&mate, &record), 150);

        // --> <--
        // This is the example given in _Sequence Alignment/Map Format Specification_ (2021-06-03)
        // § 1.4.9 "TLEN" (footnote 14).
        let record = Record::builder()
            .set_alignment_start(Position::try_from(100)?)
            .set_read_length(50)
            .build();

        let mate = Record::builder()
            .set_bam_flags(Flags::REVERSE_COMPLEMENTED)
            .set_alignment_start(Position::try_from(200)?)
            .set_read_length(50)
            .build();

        assert_eq!(calculate_template_size(&record, &mate), 150);
        assert_eq!(calculate_template_size(&mate, &record), 150);

        // <-- -->
        let record = Record::builder()
            .set_bam_flags(Flags::REVERSE_COMPLEMENTED)
            .set_alignment_start(Position::try_from(100)?)
            .set_read_length(50)
            .build();

        let mate = Record::builder()
            .set_alignment_start(Position::try_from(200)?)
            .set_read_length(50)
            .build();

        assert_eq!(calculate_template_size(&record, &mate), 150);
        assert_eq!(calculate_template_size(&mate, &record), 150);

        // <-- <--
        let record = Record::builder()
            .set_bam_flags(Flags::REVERSE_COMPLEMENTED)
            .set_alignment_start(Position::try_from(100)?)
            .set_read_length(50)
            .build();

        let mate = Record::builder()
            .set_bam_flags(Flags::REVERSE_COMPLEMENTED)
            .set_alignment_start(Position::try_from(200)?)
            .set_read_length(50)
            .build();

        assert_eq!(calculate_template_size(&record, &mate), 150);
        assert_eq!(calculate_template_size(&mate, &record), 150);

        Ok(())
    }

    #[test]
    fn test_resolve_bases() -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZeroUsize;

        use sam::{
            alignment::record_buf::Sequence,
            header::record::value::map::{self, Map},
        };

        use crate::{
            container::block::ContentType,
            record::{Feature, Features},
        };

        const SQ0_LENGTH: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        let start = Position::try_from(1)?;
        let end = Position::try_from(2)?;
        let sequence = fasta::record::Sequence::from(b"ACGT".to_vec());
        let reference_md5 = builder::calculate_normalized_sequence_digest(&sequence[start..=end]);

        let reference_sequence_repository = fasta::Repository::new(vec![fasta::Record::new(
            fasta::record::Definition::new("sq0", None),
            sequence,
        )]);

        let header = sam::Header::builder()
            .add_reference_sequence(
                "sq0".parse()?,
                Map::<map::ReferenceSequence>::new(SQ0_LENGTH),
            )
            .build();

        let compression_header = CompressionHeader::builder().build();

        let slice = Slice {
            header: Header::builder()
                .set_reference_sequence_context(ReferenceSequenceContext::some(0, start, end))
                .set_reference_md5(reference_md5)
                .build(),
            core_data_block: Block::builder()
                .set_content_type(ContentType::CoreData)
                .build(),
            external_blocks: vec![Block::builder()
                .set_content_type(ContentType::ExternalData)
                .build()],
        };

        let mut records = [Record::builder()
            .set_id(1)
            .set_bam_flags(sam::alignment::record_buf::Flags::default())
            .set_reference_sequence_id(0)
            .set_read_length(2)
            .set_alignment_start(Position::MIN)
            .set_features(Features::from(vec![Feature::Bases(
                Position::MIN,
                vec![b'A', b'C'],
            )]))
            .build()];

        resolve_bases(
            &reference_sequence_repository,
            &header,
            &compression_header,
            &slice,
            &mut records,
        )?;

        let actual: Vec<_> = records.into_iter().map(|r| r.bases).collect();
        let expected = [Sequence::from(vec![b'A', b'C'])];
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_resolve_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        use sam::alignment::record_buf::QualityScores;

        use crate::record::{Feature, Features};

        let mut records = [
            Record::builder()
                .set_id(1)
                .set_bam_flags(sam::alignment::record_buf::Flags::empty())
                .set_read_length(2)
                .set_features(Features::from(vec![Feature::Scores(
                    Position::try_from(1)?,
                    vec![8, 13],
                )]))
                .build(),
            Record::builder().set_id(2).build(),
            Record::builder()
                .set_id(3)
                .set_flags(Flags::QUALITY_SCORES_STORED_AS_ARRAY)
                .set_read_length(2)
                .set_quality_scores(QualityScores::from(vec![21, 34]))
                .build(),
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
