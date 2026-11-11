mod header;
pub mod records;

use std::{borrow::Cow, io, iter, sync::Arc};

use noodles_core::Position;
use noodles_fasta as fasta;
use noodles_sam::{self as sam, alignment::Record as _};

use self::{
    header::read_header,
    records::{ExternalDataReaders, Records},
};
use super::read_block_as;
use crate::{
    Record, calculate_normalized_sequence_digest,
    container::{
        CompressionHeader, ReferenceSequenceContext,
        block::{self, ContentType},
        slice::Header,
    },
    io::BitReader,
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
    pub fn decode_blocks(
        &self,
    ) -> io::Result<(Cow<'c, [u8]>, Vec<(block::ContentId, Cow<'c, [u8]>)>)> {
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
        external_data_srcs: &'c [(block::ContentId, Cow<'c, [u8]>)],
    ) -> io::Result<Vec<Record<'c>>> {
        let core_data_reader = BitReader::new(core_data_src);

        let mut external_data_readers = ExternalDataReaders::new();

        for (block_content_id, src) in external_data_srcs {
            external_data_readers.insert(*block_content_id, src);
        }

        let reference_sequence_context = self.header.reference_sequence_context();
        let initial_id = self.header.record_counter();

        let mut reader = Records::new(
            compression_header,
            core_data_reader,
            external_data_readers,
            reference_sequence_context,
            initial_id,
        )?;

        let slice_reference_sequence = get_slice_reference_sequence(
            &reference_sequence_repository.clone(),
            header,
            compression_header,
            &self.header,
            external_data_srcs,
        )?;

        let substitution_matrix = compression_header.preservation_map().substitution_matrix();

        let mut records = vec![Record::default(); self.header.record_count()];

        for record in &mut records {
            reader.read_record(record)?;

            record.header = Some(header);

            if !record.bam_flags.is_unmapped() && !record.cram_flags.sequence_is_missing() {
                record.reference_sequence = if reference_sequence_context.is_many() {
                    get_record_reference_sequence(&reference_sequence_repository, header, record)?
                } else {
                    slice_reference_sequence.clone()
                };

                record.substitution_matrix = substitution_matrix.clone();
            }
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
                mate.name = record.name.clone();
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
        let segments = iter::successors(Some(i), |&j| mate_indices[j]).map(|j| &records[j]);
        let span = calculate_template_span(segments);

        let mut next_index = Some(i);

        while let Some(j) = next_index {
            records[j].template_length = span.map_or(0, |s| s.signed(&records[j]));
            next_index = mate_indices[j].take();
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

fn set_mate_chunk(
    record_bam_flags: &mut sam::alignment::record::Flags,
    record_mate_reference_sequence_id: &mut Option<usize>,
    record_mate_alignment_start: &mut Option<Position>,
    mate_bam_flags: sam::alignment::record::Flags,
    mate_reference_sequence_id: Option<usize>,
    mate_alignment_start: Option<Position>,
) {
    if mate_bam_flags.is_reverse_complemented() {
        *record_bam_flags |= sam::alignment::record::Flags::MATE_REVERSE_COMPLEMENTED;
    }

    if mate_bam_flags.is_unmapped() {
        *record_bam_flags |= sam::alignment::record::Flags::MATE_UNMAPPED;
    }

    *record_mate_reference_sequence_id = mate_reference_sequence_id;
    *record_mate_alignment_start = mate_alignment_start;
}

#[derive(Clone, Copy, Debug)]
struct TemplateSpan {
    min_start: Position,
    min_start_count: usize,
    max_end: Position,
    max_end_count: usize,
    len: i32,
}

impl TemplateSpan {
    fn signed(&self, record: &Record<'_>) -> i32 {
        let is_positive = match (record.alignment_start, record.alignment_end()) {
            (Some(start), Some(end)) if start == self.min_start => {
                if self.min_start_count == 1 || end < self.max_end {
                    true
                } else if self.max_end_count > 1 {
                    record.bam_flags.is_first_segment()
                } else {
                    false
                }
            }
            _ => false,
        };

        if is_positive { self.len } else { -self.len }
    }
}

// _Sequence Alignment/Map Format Specification_ (2025-08-12) § 1.4.9 "TLEN"
fn calculate_template_span<'r, 'c: 'r, I>(records: I) -> Option<TemplateSpan>
where
    I: IntoIterator<Item = &'r Record<'c>>,
{
    let mut iter = records.into_iter();

    let first_record = iter.next()?;

    if first_record.bam_flags.is_unmapped() {
        return None;
    }

    let reference_sequence_id = first_record.reference_sequence_id?;

    let mut min_start = first_record.alignment_start?;
    let mut min_start_count = 1;

    let mut max_end = first_record.alignment_end()?;
    let mut max_end_count = 1;

    for record in iter {
        if record.bam_flags.is_unmapped()
            || record.reference_sequence_id != Some(reference_sequence_id)
        {
            return None;
        }

        let start = record.alignment_start?;
        let end = record.alignment_end()?;

        if start < min_start {
            min_start = start;
            min_start_count = 1;
        } else if start == min_start {
            min_start_count += 1;
        }

        if end > max_end {
            max_end = end;
            max_end_count = 1;
        } else if end == max_end {
            max_end_count += 1;
        }
    }

    // "...the absolute value of TLEN equals the distance between the mapped end of the template
    // and the mapped start of the template, inclusively..."
    let n = usize::from(min_start).abs_diff(usize::from(max_end)) + 1;
    let len = i32::try_from(n).expect("invalid template length");

    Some(TemplateSpan {
        min_start,
        min_start_count,
        max_end,
        max_end_count,
        len,
    })
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum ReferenceSequence<'c> {
    Embedded {
        reference_start: Position,
        sequence: &'c [u8],
    },
    External {
        sequence: Arc<fasta::record::Sequence>,
    },
}

fn get_slice_reference_sequence<'c>(
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    compression_header: &CompressionHeader,
    slice_header: &Header,
    external_data_srcs: &'c [(block::ContentId, Cow<'c, [u8]>)],
) -> io::Result<Option<ReferenceSequence<'c>>> {
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

        // § 8.5 "Slice header block" (2024-09-04): "MD5sums should not be validated if the stored
        // checksum is all-zero."
        if let Some(expected_md5) = slice_header.reference_md5() {
            let interval = context.alignment_start()..=context.alignment_end();
            let subsequence = &sequence[interval];
            validate_sequence(subsequence, expected_md5)?;
        }

        Ok(Some(ReferenceSequence::External { sequence }))
    } else if let Some(block_content_id) = embedded_reference_bases_block_content_id {
        let sequence = external_data_srcs
            .iter()
            .find(|(id, _)| *id == block_content_id)
            .map(|(_, src)| src)
            .expect("invalid block content ID");

        Ok(Some(ReferenceSequence::Embedded {
            reference_start: context.alignment_start(),
            sequence,
        }))
    } else {
        Ok(None)
    }
}

fn get_record_reference_sequence<'c>(
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    record: &Record<'c>,
) -> io::Result<Option<ReferenceSequence<'c>>> {
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

fn validate_sequence(sequence: &[u8], expected_checksum: &[u8; 16]) -> io::Result<()> {
    let actual_checksum = calculate_normalized_sequence_digest(sequence);

    if &actual_checksum == expected_checksum {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "reference sequence checksum mismatch: expected {expected_checksum:?}, got {actual_checksum:?}"
            ),
        ))
    }
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
                bam_flags: sam::alignment::record::Flags::SEGMENTED,
                cram_flags: Flags::MATE_IS_DOWNSTREAM,
                reference_sequence_id: Some(2),
                read_length: 4,
                alignment_start: Position::new(5),
                mate_distance: Some(0),
                ..Default::default()
            },
            Record {
                id: 2,
                bam_flags: sam::alignment::record::Flags::SEGMENTED,
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
                bam_flags: sam::alignment::record::Flags::SEGMENTED,
                reference_sequence_id: Some(2),
                read_length: 4,
                alignment_start: Position::new(13),
                ..Default::default()
            },
            Record {
                id: 5,
                bam_flags: sam::alignment::record::Flags::SEGMENTED,
                cram_flags: Flags::MATE_IS_DOWNSTREAM,
                reference_sequence_id: Some(2),
                read_length: 4,
                alignment_start: Position::new(11),
                mate_distance: Some(0),
                ..Default::default()
            },
            Record {
                id: 6,
                bam_flags: sam::alignment::record::Flags::SEGMENTED,
                cram_flags: Flags::MATE_IS_DOWNSTREAM,
                reference_sequence_id: Some(2),
                read_length: 4,
                alignment_start: Position::new(5),
                mate_distance: Some(0),
                ..Default::default()
            },
            Record {
                id: 7,
                bam_flags: sam::alignment::record::Flags::SEGMENTED,
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

        let name_5 = b"5".as_bstr();

        assert_eq!(records[4].name(), Some(name_5));
        assert_eq!(
            records[4].mate_reference_sequence_id,
            records[5].reference_sequence_id
        );
        assert_eq!(records[4].mate_alignment_start, records[5].alignment_start);
        assert_eq!(records[4].template_length, -12);

        assert_eq!(records[5].name(), Some(name_5));
        assert_eq!(
            records[5].mate_reference_sequence_id,
            records[6].reference_sequence_id
        );
        assert_eq!(records[5].mate_alignment_start, records[6].alignment_start);
        assert_eq!(records[5].template_length, 12);

        assert_eq!(records[6].name(), Some(name_5));
        assert_eq!(
            records[6].mate_reference_sequence_id,
            records[4].reference_sequence_id
        );
        assert_eq!(records[6].mate_alignment_start, records[4].alignment_start);
        assert_eq!(records[6].template_length, -12);

        Ok(())
    }

    #[test]
    fn test_calculate_template_span() {
        use sam::alignment::record::Flags;

        // --> -->
        let record = Record {
            bam_flags: Flags::SEGMENTED,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            bam_flags: Flags::SEGMENTED,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(200),
            read_length: 50,
            ..Default::default()
        };

        assert_eq!(
            calculate_template_span([&record, &mate])
                .map(|span| (span.signed(&record), span.signed(&mate))),
            Some((150, -150))
        );

        // --> <--
        // This is the example given in _Sequence Alignment/Map Format Specification_ (2025-08-12)
        // § 1.4.9 "TLEN" (footnote 15).
        let record = Record {
            bam_flags: Flags::SEGMENTED,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            bam_flags: Flags::SEGMENTED | Flags::REVERSE_COMPLEMENTED,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(200),
            read_length: 50,
            ..Default::default()
        };

        assert_eq!(
            calculate_template_span([&record, &mate])
                .map(|span| (span.signed(&record), span.signed(&mate))),
            Some((150, -150))
        );

        // <-- -->
        let record = Record {
            bam_flags: Flags::SEGMENTED | Flags::REVERSE_COMPLEMENTED,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            bam_flags: Flags::SEGMENTED,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(200),
            read_length: 50,
            ..Default::default()
        };

        assert_eq!(
            calculate_template_span([&record, &mate])
                .map(|span| (span.signed(&record), span.signed(&mate))),
            Some((150, -150))
        );

        // <-- <--
        let record = Record {
            bam_flags: Flags::SEGMENTED | Flags::REVERSE_COMPLEMENTED,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            bam_flags: Flags::SEGMENTED | Flags::REVERSE_COMPLEMENTED,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(200),
            read_length: 50,
            ..Default::default()
        };

        assert_eq!(
            calculate_template_span([&record, &mate])
                .map(|span| (span.signed(&record), span.signed(&mate))),
            Some((150, -150))
        );

        // same starts, different ends
        let record = Record {
            bam_flags: Flags::SEGMENTED | Flags::FIRST_SEGMENT,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            bam_flags: Flags::SEGMENTED | Flags::LAST_SEGMENT,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(100),
            read_length: 100,
            ..Default::default()
        };

        assert_eq!(
            calculate_template_span([&record, &mate])
                .map(|span| (span.signed(&record), span.signed(&mate))),
            Some((100, -100))
        );

        // record == mate
        let record = Record {
            bam_flags: Flags::SEGMENTED | Flags::FIRST_SEGMENT,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            bam_flags: Flags::SEGMENTED | Flags::LAST_SEGMENT,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        assert_eq!(
            calculate_template_span([&record, &mate])
                .map(|span| (span.signed(&record), span.signed(&mate))),
            Some((50, -50))
        );

        // placed unmapped mate
        let record = Record {
            bam_flags: Flags::SEGMENTED,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            bam_flags: Flags::SEGMENTED | Flags::UNMAPPED,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(200),
            read_length: 50,
            ..Default::default()
        };

        assert!(calculate_template_span([&record, &mate]).is_none());

        // reference sequence ID mismatch
        let record = Record {
            bam_flags: Flags::SEGMENTED,
            reference_sequence_id: Some(0),
            alignment_start: Position::new(100),
            read_length: 50,
            ..Default::default()
        };

        let mate = Record {
            bam_flags: Flags::SEGMENTED,
            reference_sequence_id: Some(1),
            alignment_start: Position::new(200),
            read_length: 50,
            ..Default::default()
        };

        assert!(calculate_template_span([&record, &mate]).is_none());

        // no alignment start position
        let record = Record::default();
        assert!(calculate_template_span([&record, &record]).is_none());

        // no segments
        assert!(calculate_template_span(iter::empty()).is_none());
    }
}
