mod header;
pub mod records;

use std::{borrow::Cow, io};

use noodles_core::Position;
use noodles_fasta as fasta;
use noodles_sam::{self as sam, alignment::Record as _};

use self::{
    header::read_header,
    records::{ExternalDataReaders, Records},
};
use super::read_block_as;
use crate::{
    calculate_normalized_sequence_digest,
    container::{
        block::{self, ContentType},
        slice::Header,
        CompressionHeader, ReferenceSequenceContext,
    },
    io::BitReader,
    record::Feature,
    Record,
};

/// A container slice.
///
/// A slice contains a header, a core data block, and one or more external blocks. This is where
/// the CRAM records are stored.
pub struct Slice<'c> {
    header: Header,
    src: &'c [u8],
}

impl<'c> Slice<'c> {
    pub(crate) fn header(&self) -> &Header {
        &self.header
    }

    #[allow(clippy::type_complexity)]
    pub fn decode_blocks(&self) -> io::Result<(Vec<u8>, Vec<(block::ContentId, Vec<u8>)>)> {
        let mut src = self.src;

        let block = read_block_as(&mut src, ContentType::CoreData)?;
        let core_data_src = block.decode()?;

        let external_data_block_count = self.header.block_count() - 1;
        let external_data_srcs = (0..external_data_block_count)
            .map(|_| {
                let block = read_block_as(&mut src, ContentType::ExternalData)?;
                block.decode().map(|src| (block.content_id, src))
            })
            .collect::<io::Result<_>>()?;

        Ok((core_data_src, external_data_srcs))
    }

    /// Reads and returns a list of raw records in this slice.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_cram::{self as cram, io::reader::Container};
    /// use noodles_fasta as fasta;
    ///
    /// let data = [];
    /// let mut reader = cram::io::Reader::new(&data[..]);
    /// let header = reader.read_header()?;
    ///
    /// let mut container = Container::default();
    ///
    /// while reader.read_container(&mut container)? != 0 {
    ///     let compression_header = container.compression_header()?;
    ///
    ///     for result in container.slices() {
    ///         let slice = result?;
    ///
    ///         let (core_data_src, external_data_srcs) = slice.decode_blocks()?;
    ///
    ///         let records = slice.records(
    ///             fasta::Repository::default(),
    ///             &header,
    ///             &compression_header,
    ///             &core_data_src,
    ///             &external_data_srcs,
    ///         )?;
    ///
    ///         // ...
    ///     }
    /// }
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn records<'h: 'c, 'ch: 'c>(
        &self,
        reference_sequence_repository: fasta::Repository,
        header: &'h sam::Header,
        compression_header: &'ch CompressionHeader,
        core_data_src: &'c [u8],
        external_data_srcs: &'c [(block::ContentId, Vec<u8>)],
    ) -> io::Result<Vec<Record>> {
        let core_data_reader = BitReader::new(core_data_src);

        let mut external_data_readers = ExternalDataReaders::new();

        for (block_content_id, src) in external_data_srcs {
            external_data_readers.insert(*block_content_id, src);
        }

        let reference_sequence_context = self.header.reference_sequence_context();

        let mut reader = Records::new(
            compression_header,
            core_data_reader,
            external_data_readers,
            reference_sequence_context,
        );

        let slice_reference_sequence = get_slice_reference_sequence(
            &reference_sequence_repository.clone(),
            header,
            compression_header,
            &self.header,
            external_data_srcs,
        )?;

        let record_count = self.header.record_count();

        let mut records = vec![Record::default(); record_count];

        let start_id = self.header.record_counter();
        let end_id = start_id + (record_count as u64);
        let ids = start_id..end_id;

        let substitution_matrix = compression_header.preservation_map().substitution_matrix();

        for (id, record) in ids.zip(&mut records) {
            record.id = id;
            reader.read_record(record)?;

            record.header = Some(header);

            record.reference_sequence = if reference_sequence_context.is_many() {
                get_record_reference_sequence(&reference_sequence_repository, header, record)?
            } else {
                slice_reference_sequence.clone()
            };

            record.substitution_matrix = substitution_matrix.clone();
        }

        resolve_mates(&mut records)?;

        Ok(records)
    }
}

pub fn read_slice<'c>(src: &mut &'c [u8]) -> io::Result<Slice<'c>> {
    let header = read_header(src)?;
    Ok(Slice { header, src })
}

fn resolve_mates(records: &mut [Record]) -> io::Result<()> {
    let mut mate_indices: Vec<_> = records
        .iter()
        .enumerate()
        .map(|(i, record)| record.mate_distance.map(|len| i + len + 1))
        .collect();

    for i in 0..records.len() {
        let record = &mut records[i];

        if record.name.is_none() {
            let name = record.id.to_string().into_bytes();
            record.name = Some(Cow::from(name));
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

            if mate.name.is_none() {
                mate.name = record.name.as_ref().map(|name| name.to_vec().into());
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
        let template_length = calculate_template_length(record, mate);
        records[i].template_length = template_length;

        let mut j = i;

        while let Some(mate_index) = mate_indices[j] {
            let record = &mut records[mate_index];
            record.template_length = -template_length;
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
        mate.bam_flags,
        mate.reference_sequence_id,
        mate.alignment_start,
    );
}

fn calculate_template_length(record: &Record, mate: &Record) -> i32 {
    calculate_template_length_chunk(
        record.alignment_start,
        record.read_length,
        &record.features,
        mate.alignment_start,
        mate.read_length,
        &mate.features,
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
fn calculate_template_length_chunk(
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
    let len = if start > end {
        start - end + 1
    } else {
        end - start + 1
    };

    i32::try_from(len).expect("invalid tempalte length")
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum ReferenceSequence {
    Embedded {
        reference_start: Position,
        sequence: fasta::record::Sequence,
    },
    External {
        sequence: fasta::record::Sequence,
    },
}

fn get_slice_reference_sequence(
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    compression_header: &CompressionHeader,
    slice_header: &Header,
    external_data_srcs: &[(block::ContentId, Vec<u8>)],
) -> io::Result<Option<ReferenceSequence>> {
    let reference_sequence_context = slice_header.reference_sequence_context();

    let ReferenceSequenceContext::Some(context) = reference_sequence_context else {
        return Ok(None);
    };

    let external_reference_sequence_is_required = compression_header
        .preservation_map()
        .external_reference_sequence_is_required();

    let embedded_reference_bases_block_content_id =
        slice_header.embedded_reference_bases_block_content_id();

    if external_reference_sequence_is_required {
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

        Ok(Some(ReferenceSequence::External { sequence }))
    } else if let Some(block_content_id) = embedded_reference_bases_block_content_id {
        let src = external_data_srcs
            .iter()
            .find(|(id, _)| *id == block_content_id)
            .map(|(_, src)| src)
            .expect("invalid block content ID");

        Ok(Some(ReferenceSequence::Embedded {
            reference_start: context.alignment_start(),
            sequence: fasta::record::Sequence::from(src.clone()),
        }))
    } else {
        Ok(None)
    }
}

fn get_record_reference_sequence(
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    record: &Record<'_>,
) -> io::Result<Option<ReferenceSequence>> {
    if record.bam_flags.is_unmapped() {
        return Ok(None);
    }

    let reference_sequence_name = record
        .reference_sequence(header)
        .transpose()?
        .map(|(name, _)| name)
        .expect("invalid reference sequence ID");

    let sequence = reference_sequence_repository
        .get(reference_sequence_name)
        .transpose()?
        .expect("invalid reference sequence name");

    Ok(Some(ReferenceSequence::External { sequence }))
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;

    use super::*;
    use crate::record::Flags;

    #[test]
    fn test_resolve_mates() -> io::Result<()> {
        let mut records = vec![
            Record {
                id: 1,
                cram_flags: Flags::MATE_IS_DOWNSTREAM,
                reference_sequence_id: Some(2),
                read_length: 4,
                alignment_start: Position::new(5),
                mate_distance: Some(0),
                ..Default::default()
            },
            Record {
                id: 2,
                cram_flags: Flags::MATE_IS_DOWNSTREAM,
                reference_sequence_id: Some(2),
                read_length: 4,
                alignment_start: Position::new(8),
                mate_distance: Some(1),
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
            records[0].mate_reference_sequence_id,
            records[1].reference_sequence_id
        );
        assert_eq!(records[0].mate_alignment_start, records[1].alignment_start);
        assert_eq!(records[0].template_length, 12);

        assert_eq!(records[1].name(), Some(name_1));
        assert_eq!(
            records[1].mate_reference_sequence_id,
            records[3].reference_sequence_id
        );
        assert_eq!(records[1].mate_alignment_start, records[3].alignment_start);
        assert_eq!(records[1].template_length, -12);

        let name_3 = b"3".as_bstr();
        assert_eq!(records[2].name(), Some(name_3));

        assert_eq!(records[3].name(), Some(name_1));
        assert_eq!(
            records[3].mate_reference_sequence_id,
            records[0].reference_sequence_id
        );
        assert_eq!(records[3].mate_alignment_start, records[0].alignment_start);
        assert_eq!(records[3].template_length, -12);

        Ok(())
    }

    #[test]
    fn test_calculate_template_length() {
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

        assert_eq!(calculate_template_length(&record, &mate), 150);
        assert_eq!(calculate_template_length(&mate, &record), 150);

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

        assert_eq!(calculate_template_length(&record, &mate), 150);
        assert_eq!(calculate_template_length(&mate, &record), 150);

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

        assert_eq!(calculate_template_length(&record, &mate), 150);
        assert_eq!(calculate_template_length(&mate, &record), 150);

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

        assert_eq!(calculate_template_length(&record, &mate), 150);
        assert_eq!(calculate_template_length(&mate, &record), 150);

        // No alignment start position.
        let record = Record::default();
        assert_eq!(calculate_template_length(&record, &record), 0);
    }
}
