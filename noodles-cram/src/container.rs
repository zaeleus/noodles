pub mod block;
pub mod compression_header;
mod header;
pub mod reference_sequence_id;
pub mod slice;

pub use self::{
    block::Block, compression_header::CompressionHeader, header::Header,
    reference_sequence_id::ReferenceSequenceId, slice::Slice,
};

use std::{cmp, io};

use super::{
    num::Itf8, writer, writer::compression_header::write_compression_header, DataContainer,
};

#[derive(Debug, Default)]
pub struct Container {
    header: Header,
    blocks: Vec<Block>,
}

impl Container {
    /// Creates an EOF container.
    pub fn eof() -> Self {
        Self::new(Header::eof(), vec![Block::eof()])
    }

    pub fn try_from_data_container(data_container: &DataContainer) -> io::Result<Self> {
        let mut buf = Vec::new();
        write_compression_header(&mut buf, data_container.compression_header())?;

        let block = Block::new(
            block::CompressionMethod::None,
            block::ContentType::CompressionHeader,
            0, // FIXME
            buf.len() as Itf8,
            buf,
            0,
        );

        let mut blocks = vec![block];
        let mut landmarks = Vec::new();

        let mut container_reference_sequence_id: Option<ReferenceSequenceId> = None;

        let mut container_alignment_start = i32::MAX;
        let mut container_alignment_end = 1;

        let mut container_record_count = 0;

        for slice in data_container.slices() {
            let slice_header = slice.header();

            if let Some(reference_sequence_id) = container_reference_sequence_id {
                if !reference_sequence_id.is_many()
                    && reference_sequence_id != slice.header().reference_sequence_id()
                {
                    container_reference_sequence_id = Some(ReferenceSequenceId::Many);
                }
            } else {
                container_reference_sequence_id = Some(slice.header().reference_sequence_id());
            }

            container_alignment_start =
                cmp::min(container_alignment_start, slice_header.alignment_start());

            let slice_alignment_end =
                slice_header.alignment_start() + slice_header.alignment_span() - 1;
            container_alignment_end = cmp::max(container_alignment_end, slice_alignment_end);

            container_record_count += slice_header.record_count();

            let mut slice_len = 0;

            let mut slice_header_buf = Vec::new();
            writer::slice::write_header(&mut slice_header_buf, slice.header())?;

            let slice_header_block = Block::new(
                block::CompressionMethod::None,
                block::ContentType::SliceHeader,
                0, // FIXME,
                slice_header_buf.len() as Itf8,
                slice_header_buf,
                0,
            );

            slice_len += slice_header_block.len() as Itf8;
            blocks.push(slice_header_block);

            blocks.push(slice.core_data_block().clone());
            slice_len += slice.core_data_block().len() as Itf8;

            for external_block in slice.external_blocks() {
                blocks.push(external_block.clone());
                slice_len += external_block.len() as Itf8;
            }

            let last_landmark = landmarks.last().unwrap_or(&0);
            let landmark = last_landmark + slice_len;
            landmarks.push(landmark);
        }

        let len = blocks.iter().map(|b| b.len() as i32).sum();

        let container_alignment_span = container_alignment_end - container_alignment_start + 1;

        let header = Header::builder()
            .set_length(len)
            .set_reference_sequence_id(
                container_reference_sequence_id.expect("no slices in builder"),
            )
            .set_start_position(container_alignment_start)
            .set_alignment_span(container_alignment_span)
            .set_record_count(container_record_count)
            // TODO: set record counter
            // TODO: set base count
            .set_block_count(blocks.len() as Itf8)
            .set_landmarks(landmarks)
            .build();

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

    pub fn is_eof(&self) -> bool {
        self.header.is_eof()
    }
}
