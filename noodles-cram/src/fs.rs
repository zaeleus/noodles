//! CRAM filesystem operations.

use std::{cmp, collections::HashMap, fs::File, io, path::Path};

use noodles_core::Position;

use super::{
    container::{slice, CompressionHeader},
    crai,
    io::{reader::container::Slice, Reader},
};
use crate::io::reader::Container;

/// Indexes a CRAM file.
///
/// # Examples
///
/// ```no_run
/// use noodles_cram as cram;
/// let index = cram::fs::index("sample.cram")?;
/// # Ok::<(), std::io::Error>(())
/// ```
pub fn index<P>(src: P) -> io::Result<crai::Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    reader.read_header()?;

    let mut index = Vec::new();

    let mut container = Container::default();
    let mut container_position = reader.position()?;

    loop {
        let container_len = match reader.read_container(&mut container)? {
            0 => break,
            n => n,
        };

        let compression_header = container.compression_header()?;

        let landmarks = container.header().landmarks();
        let slice_count = landmarks.len();

        for (i, result) in container.slices().enumerate() {
            let slice = result?;
            let landmark = landmarks[i];

            let slice_length = if i < slice_count - 1 {
                landmarks[i + 1] - landmark
            } else {
                container_len - landmark
            };

            push_index_records(
                &mut index,
                &compression_header,
                &slice,
                container_position,
                landmark as u64,
                slice_length as u64,
            )?;
        }

        container_position = reader.position()?;
    }

    Ok(index)
}

fn push_index_records(
    index: &mut crai::Index,
    compression_header: &CompressionHeader,
    slice: &Slice,
    container_position: u64,
    landmark: u64,
    slice_length: u64,
) -> io::Result<()> {
    if slice.header().reference_sequence_context().is_many() {
        push_index_records_for_multi_reference_slice(
            index,
            compression_header,
            slice,
            container_position,
            landmark,
            slice_length,
        )
    } else {
        push_index_record_for_single_reference_slice(
            index,
            slice.header(),
            container_position,
            landmark,
            slice_length,
        )
    }
}

#[derive(Debug)]
struct SliceReferenceSequenceAlignmentRangeInclusive {
    start: Option<Position>,
    end: Option<Position>,
}

impl Default for SliceReferenceSequenceAlignmentRangeInclusive {
    fn default() -> Self {
        Self {
            start: Position::new(usize::MAX),
            end: None,
        }
    }
}

fn push_index_records_for_multi_reference_slice(
    index: &mut crai::Index,
    compression_header: &CompressionHeader,
    slice: &Slice,
    container_position: u64,
    landmark: u64,
    slice_length: u64,
) -> io::Result<()> {
    let mut reference_sequence_ids: HashMap<
        Option<usize>,
        SliceReferenceSequenceAlignmentRangeInclusive,
    > = HashMap::new();

    let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

    for record in slice.records(compression_header, &core_data_src, &external_data_srcs)? {
        let reference_sequence_id = record.reference_sequence_id();

        let range = reference_sequence_ids
            .entry(reference_sequence_id)
            .or_default();

        let alignment_start = record.alignment_start();
        range.start = cmp::min(range.start, alignment_start);

        let alignment_end = record.alignment_end();
        range.end = cmp::max(range.end, alignment_end);
    }

    let mut sorted_reference_sequence_ids: Vec<_> =
        reference_sequence_ids.keys().copied().collect();
    sorted_reference_sequence_ids.sort_unstable();

    for reference_sequence_id in sorted_reference_sequence_ids {
        let (alignment_start, alignment_span) = if reference_sequence_id.is_some() {
            let range = &reference_sequence_ids[&reference_sequence_id];

            if let (Some(start), Some(end)) = (range.start, range.end) {
                let span = usize::from(end) - usize::from(start) + 1;
                (Some(start), span)
            } else {
                todo!("unhandled interval: {:?}", range);
            }
        } else {
            (None, 0)
        };

        let record = crai::Record::new(
            reference_sequence_id,
            alignment_start,
            alignment_span,
            container_position,
            landmark,
            slice_length,
        );

        index.push(record);
    }

    Ok(())
}

fn push_index_record_for_single_reference_slice(
    index: &mut crai::Index,
    slice_header: &slice::Header,
    container_position: u64,
    landmark: u64,
    slice_length: u64,
) -> io::Result<()> {
    use crate::container::ReferenceSequenceContext;

    let (reference_sequence_id, alignment_start, alignment_span) =
        match slice_header.reference_sequence_context() {
            ReferenceSequenceContext::Some(context) => {
                let reference_sequence_id = Some(context.reference_sequence_id());
                let alignment_start = Some(context.alignment_start());
                let alignment_span = context.alignment_span();
                (reference_sequence_id, alignment_start, alignment_span)
            }
            ReferenceSequenceContext::None => (None, None, 0),
            ReferenceSequenceContext::Many => unreachable!(),
        };

    let record = crai::Record::new(
        reference_sequence_id,
        alignment_start,
        alignment_span,
        container_position,
        landmark,
        slice_length,
    );

    index.push(record);

    Ok(())
}
