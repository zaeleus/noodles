#![warn(missing_docs)]

//! **noodles-cram** handles the reading and writing of the CRAM format.

#[cfg(feature = "async")]
#[allow(dead_code)]
mod r#async;

mod bit_reader;
mod bit_writer;
pub(crate) mod container;
pub mod crai;
mod data_container;
pub mod file_definition;
mod huffman;
mod num;
mod rans;
pub mod reader;
pub mod record;
pub(crate) mod writer;

pub use self::{file_definition::FileDefinition, reader::Reader, record::Record, writer::Writer};

#[cfg(feature = "async")]
pub use self::r#async::Reader as AsyncReader;

pub(crate) use self::{
    bit_reader::BitReader, bit_writer::BitWriter, container::Container,
    data_container::DataContainer,
};

use std::{cmp, collections::HashMap, convert::TryFrom, fs::File, io, path::Path};

use noodles_bam as bam;

static MAGIC_NUMBER: &[u8] = b"CRAM";

/// Indexes a CRAM file.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// use noodles_cram as cram;
/// let index = cram::index("sample.cram")?;
/// # Ok::<(), io::Error>(())
/// ```
pub fn index<P>(src: P) -> io::Result<crai::Index>
where
    P: AsRef<Path>,
{
    let mut reader = File::open(src).map(Reader::new)?;
    reader.read_file_definition()?;
    reader.read_file_header()?;

    let mut index = Vec::new();
    let mut container_position = reader.position()?;

    loop {
        let container = reader.read_container()?;

        if container.is_eof() {
            break;
        }

        let container_header = container.header();
        let container_len = container_header.len();

        let landmarks = container_header.landmarks();
        let slice_count = landmarks.len();

        let data_container = DataContainer::try_from(container.clone())?;

        for (i, slice) in data_container.slices().iter().enumerate() {
            let landmark = landmarks[i];

            let slice_length = if i < slice_count - 1 {
                landmarks[i + 1] - landmark
            } else {
                container_len - landmark
            };

            push_index_records(
                &mut index,
                data_container.compression_header(),
                slice,
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
    compression_header: &data_container::CompressionHeader,
    slice: &data_container::Slice,
    container_position: u64,
    landmark: u64,
    slice_length: u64,
) -> io::Result<()> {
    if slice.header().reference_sequence_id().is_many() {
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
    start: i32,
    end: i32,
}

impl Default for SliceReferenceSequenceAlignmentRangeInclusive {
    fn default() -> Self {
        Self {
            start: i32::MAX,
            end: 1,
        }
    }
}

fn push_index_records_for_multi_reference_slice(
    index: &mut crai::Index,
    compression_header: &data_container::CompressionHeader,
    slice: &data_container::Slice,
    container_position: u64,
    landmark: u64,
    slice_length: u64,
) -> io::Result<()> {
    let mut raw_reference_sequence_ids: HashMap<
        i32,
        SliceReferenceSequenceAlignmentRangeInclusive,
    > = HashMap::new();

    for record in slice.records(compression_header)? {
        let raw_reference_sequence_id = record
            .reference_sequence_id()
            .map(i32::from)
            .unwrap_or(bam::record::reference_sequence_id::UNMAPPED);

        let range = raw_reference_sequence_ids
            .entry(raw_reference_sequence_id)
            .or_default();

        range.start = cmp::min(range.start, record.alignment_start());
        range.end = cmp::max(range.end, record.alignment_end());
    }

    let mut raw_sorted_reference_sequence_ids: Vec<_> =
        raw_reference_sequence_ids.keys().copied().collect();
    raw_sorted_reference_sequence_ids.sort_unstable();

    let reference_sequence_ids = raw_sorted_reference_sequence_ids
        .iter()
        .map(|&id| {
            let range = &raw_reference_sequence_ids[&id];

            let alignment_start = range.start;
            let alignment_span = range.end - alignment_start + 1;

            if id == bam::record::reference_sequence_id::UNMAPPED {
                Ok((None, alignment_start, alignment_span))
            } else {
                bam::record::ReferenceSequenceId::try_from(id)
                    .map(Some)
                    .map(|reference_sequence_id| {
                        (reference_sequence_id, alignment_start, alignment_span)
                    })
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            }
        })
        .collect::<Result<Vec<_>, _>>()?;

    for (reference_sequence_id, alignment_start, alignment_span) in reference_sequence_ids {
        let record = crai::Record::new(
            reference_sequence_id,
            alignment_start,
            alignment_span,
            container_position,
            landmark as u64,
            slice_length as u64,
        );

        index.push(record);
    }

    Ok(())
}

fn push_index_record_for_single_reference_slice(
    index: &mut crai::Index,
    slice_header: &data_container::slice::Header,
    container_position: u64,
    landmark: u64,
    slice_length: u64,
) -> io::Result<()> {
    let slice_reference_sequence_id = slice_header.reference_sequence_id();

    let reference_sequence_id = if slice_reference_sequence_id.is_none() {
        None
    } else {
        bam::record::ReferenceSequenceId::try_from(i32::from(slice_reference_sequence_id))
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
    };

    let record = crai::Record::new(
        reference_sequence_id,
        slice_header.alignment_start(),
        slice_header.alignment_span(),
        container_position,
        landmark,
        slice_length,
    );

    index.push(record);

    Ok(())
}
