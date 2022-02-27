pub(crate) mod builder;
pub(crate) mod header;

pub use self::{builder::Builder, header::Header};

use std::io::{self, Cursor};

use noodles_fasta as fasta;
use noodles_sam as sam;

use super::CompressionHeader;
use crate::{
    container::Block,
    record::resolve::{resolve_bases, resolve_quality_scores},
    BitReader, Record,
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
    /// let mut reader = cram::Reader::new(&data[..]);
    /// reader.read_file_definition()?;
    /// reader.read_file_header()?;
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
        use crate::reader::record::ExternalDataReaders;

        let core_data_reader = self
            .core_data_block
            .decompressed_data()
            .map(Cursor::new)
            .map(BitReader::new)?;

        let mut external_data_readers = ExternalDataReaders::new();

        for block in self.external_blocks() {
            let reader = block.decompressed_data().map(Cursor::new)?;
            external_data_readers.insert(block.content_id(), reader);
        }

        let mut record_reader = crate::reader::record::Reader::new(
            compression_header,
            core_data_reader,
            external_data_readers,
            self.header.reference_sequence_id(),
            self.header.alignment_start(),
        );

        let record_count = self.header().record_count();
        let mut records = Vec::with_capacity(record_count);

        let start_id = self.header().record_counter();
        let end_id = start_id + (record_count as i64);

        for id in start_id..end_id {
            let mut record = record_reader.read_record()?;
            record.id = id;
            records.push(record);
        }

        Ok(records)
    }

    /// Resolves records.
    ///
    /// This resolves mates, read names, bases, and quality scores.
    pub fn resolve_records(
        &self,
        reference_sequences: &[fasta::Record],
        compression_header: &CompressionHeader,
        records: Vec<Record>,
    ) -> io::Result<Vec<Record>> {
        let mut records = self.resolve_mates(records)?;
        self.resolve_bases(reference_sequences, compression_header, &mut records)?;
        self.resolve_quality_scores(&mut records);
        Ok(records)
    }

    /// Resolves mate records.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_cram as cram;
    ///
    /// let data = [];
    /// let mut reader = cram::Reader::new(&data[..]);
    /// reader.read_file_definition()?;
    /// reader.read_file_header()?;
    ///
    /// while let Some(container) = reader.read_data_container()? {
    ///     for slice in container.slices() {
    ///         let records = slice.records(container.compression_header())?;
    ///         let records = slice.resolve_mates(records)?;
    ///         // ...
    ///     }
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn resolve_mates(&self, records: Vec<Record>) -> io::Result<Vec<Record>> {
        resolve_mates(records)
    }

    fn resolve_bases(
        &self,
        reference_sequences: &[fasta::Record],
        compression_header: &CompressionHeader,
        records: &mut [Record],
    ) -> io::Result<()> {
        let embedded_reference_sequence = if let Some(block_content_id) =
            self.header().embedded_reference_bases_block_content_id()
        {
            let block = self
                .external_blocks()
                .iter()
                .find(|block| block.content_id() == block_content_id)
                .expect("invalid block content ID");

            let data = block.decompressed_data()?;
            let sequence = fasta::record::Sequence::from(Vec::from(data));

            Some(sequence)
        } else {
            None
        };

        for record in records {
            if record.bam_flags().is_unmapped() || record.flags().decode_sequence_as_unknown() {
                continue;
            }

            let reference_sequence = if compression_header
                .preservation_map()
                .is_reference_required()
            {
                let reference_sequence_id = record
                    .reference_sequence_id()
                    .map(usize::from)
                    .expect("invalid reference sequence ID");

                let record = &reference_sequences[reference_sequence_id];
                Some(record.sequence().clone())
            } else {
                embedded_reference_sequence.as_ref().cloned()
            };

            let substitution_matrix = compression_header.preservation_map().substitution_matrix();
            let alignment_start = record.alignment_start().expect("invalid alignment start");

            let bases = resolve_bases(
                reference_sequence,
                substitution_matrix,
                record.features(),
                alignment_start,
                record.read_length(),
            )?;

            record.bases = bases;
        }

        Ok(())
    }

    fn resolve_quality_scores(&self, records: &mut [Record]) {
        for record in records {
            if !record.flags().are_quality_scores_stored_as_array() {
                let quality_scores =
                    resolve_quality_scores(record.features(), record.read_length());

                record.quality_scores = quality_scores;
            }
        }
    }
}

fn resolve_mates(records: Vec<Record>) -> io::Result<Vec<Record>> {
    use std::cell::RefCell;

    let mut mate_indices = vec![None; records.len()];

    for (i, record) in records.iter().enumerate() {
        let flags = record.flags();

        if flags.has_mate_downstream() {
            let distance_to_next_fragment = record.distance_to_next_fragment() as usize;
            let mate_index = i + distance_to_next_fragment + 1;
            mate_indices[i] = Some(mate_index);
        }
    }

    let records: Vec<_> = records.into_iter().map(RefCell::new).collect();

    for (i, record_cell) in records.iter().enumerate() {
        if mate_indices[i].is_none() {
            continue;
        }

        let mut record = record_cell.borrow_mut();

        if record.read_name.is_empty() {
            let read_name = record.id().to_string().into_bytes();
            record.read_name.extend(read_name);
        }

        let mut j = i;

        while let Some(mate_index) = mate_indices[j] {
            let mut mate = records[mate_index].borrow_mut();
            set_mate(&mut record, &mut mate);
            record = mate;
            j = mate_index;
        }

        let mut mate = record_cell.borrow_mut();
        set_mate(&mut record, &mut mate);

        let template_size = calculate_template_size(&record, &mate)?;
        record.template_size = template_size;
        mate.template_size = -template_size;
    }

    Ok(records.into_iter().map(|r| r.into_inner()).collect())
}

fn set_mate(mut record: &mut Record, mate: &mut Record) {
    let mate_bam_flags = mate.bam_flags();

    if mate_bam_flags.is_reverse_complemented() {
        record.bam_bit_flags |= sam::record::Flags::MATE_REVERSE_COMPLEMENTED;
    }

    if mate_bam_flags.is_unmapped() {
        record.bam_bit_flags |= sam::record::Flags::MATE_UNMAPPED;
    }

    if mate.read_name().is_empty() {
        mate.read_name.extend(record.read_name.iter());
    }

    record.next_fragment_reference_sequence_id = mate.reference_sequence_id();
    record.next_mate_alignment_start = mate.alignment_start();
}

// _Sequence Alignment/Map Format Specification_ (2021-06-03) § 1.4.9 "TLEN"
fn calculate_template_size(record: &Record, mate: &Record) -> io::Result<i32> {
    let start = if record.bam_flags().is_reverse_complemented() {
        record
            .alignment_end()
            .transpose()?
            .map(i32::from)
            .unwrap_or_default()
    } else {
        record.alignment_start().map(i32::from).unwrap_or_default()
    };

    let end = if mate.bam_flags().is_reverse_complemented() {
        mate.alignment_end()
            .transpose()?
            .map(i32::from)
            .unwrap_or_default()
    } else {
        mate.alignment_start().map(i32::from).unwrap_or_default()
    };

    // "...the absolute value of TLEN equals the distance between the mapped end of the template
    // and the mapped start of the template, inclusively..."
    let len = (end - start).abs() + 1;

    // "The TLEN field is positive for the leftmost segment of the template, negative for the
    // rightmost, and the sign for any middle segment is undefined. If segments cover the same
    // coordinates then the choice of which is leftmost and rightmost is arbitrary..."
    if start > end {
        Ok(-len)
    } else {
        Ok(len)
    }
}

#[cfg(test)]
mod tests {
    use noodles_bam as bam;

    use super::*;

    #[test]
    fn test_resolve_mates() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::Flags;
        use bam::record::ReferenceSequenceId;

        let records = vec![
            Record::builder()
                .set_id(1)
                .set_flags(Flags::HAS_MATE_DOWNSTREAM)
                .set_reference_sequence_id(ReferenceSequenceId::try_from(2)?)
                .set_read_length(4)
                .set_alignment_start(sam::record::Position::try_from(5)?)
                .set_distance_to_next_fragment(0)
                .build(),
            Record::builder()
                .set_id(2)
                .set_flags(Flags::HAS_MATE_DOWNSTREAM)
                .set_reference_sequence_id(ReferenceSequenceId::try_from(2)?)
                .set_read_length(4)
                .set_alignment_start(sam::record::Position::try_from(8)?)
                .set_distance_to_next_fragment(1)
                .build(),
            Record::builder().set_id(3).build(),
            Record::builder()
                .set_id(4)
                .set_reference_sequence_id(ReferenceSequenceId::try_from(2)?)
                .set_read_length(4)
                .set_alignment_start(sam::record::Position::try_from(13)?)
                .build(),
        ];

        let records = resolve_mates(records)?;

        assert_eq!(records[0].read_name(), b"1");
        assert_eq!(
            records[0].next_fragment_reference_sequence_id(),
            records[1].reference_sequence_id()
        );
        assert_eq!(
            records[0].next_mate_alignment_start(),
            records[1].alignment_start(),
        );

        assert_eq!(records[1].read_name(), b"1");
        assert_eq!(
            records[1].next_fragment_reference_sequence_id(),
            records[3].reference_sequence_id()
        );
        assert_eq!(
            records[1].next_mate_alignment_start(),
            records[3].alignment_start(),
        );

        // FIXME
        // assert_eq!(records[2].read_name(), b"3");

        assert_eq!(records[3].read_name(), b"1");
        // FIXME
        /* assert_eq!(
            records[3].next_fragment_reference_sequence_id(),
            records[0].reference_sequence_id()
        );
        assert_eq!(
            records[3].next_mate_alignment_start(),
            records[0].alignment_start(),
        ); */

        Ok(())
    }

    #[test]
    fn test_calculate_template_size() -> Result<(), Box<dyn std::error::Error>> {
        use sam::record::{Flags, Position};

        // --> -->
        let record = Record::builder()
            .set_alignment_start(Position::try_from(100)?)
            .set_read_length(50)
            .build();

        let mate = Record::builder()
            .set_alignment_start(Position::try_from(200)?)
            .set_read_length(50)
            .build();

        assert_eq!(calculate_template_size(&record, &mate)?, 101);
        assert_eq!(calculate_template_size(&mate, &record)?, -101);

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

        assert_eq!(calculate_template_size(&record, &mate)?, 150);
        assert_eq!(calculate_template_size(&mate, &record)?, -150);

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

        assert_eq!(calculate_template_size(&record, &mate)?, 52);
        assert_eq!(calculate_template_size(&mate, &record)?, -52);

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

        assert_eq!(calculate_template_size(&record, &mate)?, 101);
        assert_eq!(calculate_template_size(&mate, &record)?, -101);

        Ok(())
    }
}
