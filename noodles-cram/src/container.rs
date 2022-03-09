pub mod block;
mod header;
pub mod reference_sequence_id;

pub use self::{block::Block, header::Header, reference_sequence_id::ReferenceSequenceId};

use std::{cmp, error, fmt, io};

use noodles_core::Position;
use noodles_sam as sam;

use super::{data_container::Slice, writer, DataContainer};

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Container {
    header: Header,
    blocks: Vec<Block>,
}

impl Container {
    /// Creates an EOF container.
    pub fn eof() -> Self {
        Self::new(Header::eof(), vec![Block::eof()])
    }

    pub fn try_from_data_container(
        data_container: &DataContainer,
        base_count: i64,
    ) -> io::Result<Self> {
        use super::writer::data_container::write_compression_header;

        let mut buf = Vec::new();
        write_compression_header(&mut buf, data_container.compression_header())?;

        let block = Block::builder()
            .set_content_type(block::ContentType::CompressionHeader)
            .set_uncompressed_len(buf.len())
            .set_data(buf)
            .build();

        let mut blocks = vec![block];
        let mut landmarks = Vec::new();

        let container_reference_sequence_id =
            find_container_reference_sequence_id(data_container.slices())?;

        let (container_alignment_start, container_alignment_end) =
            if container_reference_sequence_id.is_some() {
                find_container_alignment_positions(data_container.slices())?
            } else {
                (None, None)
            };

        let mut container_record_count = 0;
        let container_record_counter = data_container
            .slices()
            .first()
            .map(|s| s.header().record_counter())
            .expect("no slices in builder");

        for slice in data_container.slices() {
            let slice_header = slice.header();

            container_record_count += slice_header.record_count() as i32;

            let mut slice_len = 0;

            let mut slice_header_buf = Vec::new();
            writer::data_container::slice::write_header(&mut slice_header_buf, slice.header())?;

            let slice_header_block = Block::builder()
                .set_content_type(block::ContentType::SliceHeader)
                .set_uncompressed_len(slice_header_buf.len())
                .set_data(slice_header_buf)
                .build();

            slice_len += slice_header_block.len() as i32;
            blocks.push(slice_header_block);

            blocks.push(slice.core_data_block().clone());
            slice_len += slice.core_data_block().len() as i32;

            for external_block in slice.external_blocks() {
                blocks.push(external_block.clone());
                slice_len += external_block.len() as i32;
            }

            let last_landmark = landmarks.last().copied().unwrap_or(0);
            let landmark = last_landmark + slice_len;
            landmarks.push(landmark);
        }

        let len = blocks.iter().map(|b| b.len() as i32).sum();

        let mut builder = Header::builder()
            .set_length(len)
            .set_reference_sequence_id(container_reference_sequence_id)
            .set_record_count(container_record_count)
            .set_record_counter(container_record_counter)
            .set_base_count(base_count)
            .set_block_count(blocks.len())
            .set_landmarks(landmarks);

        if let (Some(alignment_start), Some(alignment_end)) =
            (container_alignment_start, container_alignment_end)
        {
            let alignment_span = usize::from(alignment_end) - usize::from(alignment_start) + 1;

            builder = builder
                .set_start_position(alignment_start)
                .set_alignment_span(alignment_span);
        }

        let header = builder.build();

        Ok(Self::new(header, blocks))
    }

    pub fn new(header: Header, blocks: Vec<Block>) -> Self {
        Self { header, blocks }
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn blocks(&self) -> &[Block] {
        &self.blocks
    }
}

fn find_container_reference_sequence_id(slices: &[Slice]) -> io::Result<ReferenceSequenceId> {
    assert!(!slices.is_empty());

    let first_slice = slices.first().expect("slices cannot be empty");
    let container_reference_sequence_id = first_slice.header().reference_sequence_id();

    for slice in slices.iter().skip(1) {
        let slice_reference_sequence_id = slice.header().reference_sequence_id();

        if slice_reference_sequence_id != container_reference_sequence_id {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "invalid slice reference sequence ID: expected {:?}, got {:?}",
                    container_reference_sequence_id, slice_reference_sequence_id
                ),
            ));
        }
    }

    Ok(container_reference_sequence_id)
}

fn find_container_alignment_positions(
    slices: &[Slice],
) -> io::Result<(Option<Position>, Option<Position>)> {
    assert!(!slices.is_empty());

    let mut container_alignment_start = Position::new(usize::MAX);
    let mut container_alignment_end = None;

    for slice in slices {
        let slice_header = slice.header();

        let slice_alignment_start = slice_header.alignment_start();
        container_alignment_start = cmp::min(container_alignment_start, slice_alignment_start);

        let slice_alignment_end = slice_header.alignment_end();
        container_alignment_end = cmp::max(container_alignment_end, slice_alignment_end);
    }

    Ok((container_alignment_start, container_alignment_end))
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromSamHeaderError {
    ReferenceSequenceMissingMd5Checksum,
}

impl error::Error for TryFromSamHeaderError {}

impl fmt::Display for TryFromSamHeaderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ReferenceSequenceMissingMd5Checksum => {
                f.write_str("reference sequence is missing MD5 checksum")
            }
        }
    }
}

impl TryFrom<&sam::Header> for Container {
    type Error = TryFromSamHeaderError;

    fn try_from(header: &sam::Header) -> Result<Self, Self::Error> {
        use crate::container::block::ContentType;

        validate_reference_sequences(header.reference_sequences())?;

        let header_data = header.to_string().into_bytes();
        let header_data_len = header_data.len() as i32;

        let mut data = header_data_len.to_le_bytes().to_vec();
        data.extend(header_data);

        let block = Block::builder()
            .set_content_type(ContentType::FileHeader)
            .set_uncompressed_len(data.len())
            .set_data(data)
            .build();

        let blocks = vec![block];
        let landmarks = vec![0];

        let len = blocks.iter().map(|b| b.len() as i32).sum();

        let container_header = Header::builder()
            .set_length(len)
            .set_reference_sequence_id(ReferenceSequenceId::None)
            .set_block_count(blocks.len())
            .set_landmarks(landmarks)
            .build();

        Ok(Self::new(container_header, blocks))
    }
}

fn validate_reference_sequences(
    reference_sequences: &sam::header::ReferenceSequences,
) -> Result<(), TryFromSamHeaderError> {
    for reference_sequence in reference_sequences.values() {
        if reference_sequence.md5_checksum().is_none() {
            return Err(TryFromSamHeaderError::ReferenceSequenceMissingMd5Checksum);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_sam_header_for_container() -> Result<(), Box<dyn std::error::Error>> {
        let reference_sequence = sam::header::ReferenceSequence::builder()
            .set_name("sq0".parse()?)
            .set_length(8)
            .set_md5_checksum(
                [
                    0xd7, 0xeb, 0xa3, 0x11, 0x42, 0x1b, 0xbc, 0x9d, 0x3a, 0xda, 0x44, 0x70, 0x9d,
                    0xd6, 0x15, 0x34,
                ]
                .into(),
            )
            .build()?;

        let sam_header = sam::Header::builder()
            .add_reference_sequence(reference_sequence)
            .build();

        let container = Container::try_from(&sam_header)?;

        let expected_header = Header::builder()
            .set_length(65)
            .set_reference_sequence_id(ReferenceSequenceId::None)
            .set_alignment_span(0)
            .set_record_count(0)
            .set_record_counter(0)
            .set_base_count(0)
            .set_block_count(1)
            .set_landmarks(vec![0])
            .set_crc32(0)
            .build();

        assert_eq!(container.header(), &expected_header);

        let mut data = 52_i32.to_le_bytes().to_vec();
        data.extend(sam_header.to_string().bytes());

        let expected_blocks = vec![Block::builder()
            .set_compression_method(block::CompressionMethod::None)
            .set_content_type(block::ContentType::FileHeader)
            .set_content_id(0)
            .set_uncompressed_len(56)
            .set_data(data)
            .set_crc32(0)
            .build()];

        assert_eq!(container.blocks(), expected_blocks);

        Ok(())
    }

    #[test]
    fn test_try_from_sam_header_for_container_with_missing_reference_sequence_md5_checksum(
    ) -> Result<(), Box<dyn std::error::Error>> {
        let header = sam::Header::builder()
            .add_reference_sequence(sam::header::ReferenceSequence::new("sq0".parse()?, 8)?)
            .build();

        assert_eq!(
            Container::try_from(&header),
            Err(TryFromSamHeaderError::ReferenceSequenceMissingMd5Checksum)
        );

        Ok(())
    }
}
