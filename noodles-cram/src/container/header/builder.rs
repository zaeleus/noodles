use noodles_sam as sam;

use crate::container::ReferenceSequenceId;

use super::Header;

#[derive(Debug, Default)]
pub struct Builder {
    length: i32,
    reference_sequence_id: ReferenceSequenceId,
    start_position: Option<sam::record::Position>,
    alignment_span: usize,
    record_count: i32,
    record_counter: i64,
    base_count: i64,
    block_count: usize,
    landmarks: Vec<i32>,
    crc32: u32,
}

impl Builder {
    pub fn set_length(mut self, length: i32) -> Self {
        self.length = length;
        self
    }

    pub fn set_reference_sequence_id(mut self, reference_sequence_id: ReferenceSequenceId) -> Self {
        self.reference_sequence_id = reference_sequence_id;
        self
    }

    pub fn set_start_position(mut self, start_position: sam::record::Position) -> Self {
        self.start_position = Some(start_position);
        self
    }

    pub fn set_alignment_span(mut self, alignment_span: usize) -> Self {
        self.alignment_span = alignment_span;
        self
    }

    pub fn set_record_count(mut self, record_count: i32) -> Self {
        self.record_count = record_count;
        self
    }

    pub fn set_record_counter(mut self, record_counter: i64) -> Self {
        self.record_counter = record_counter;
        self
    }

    pub fn set_base_count(mut self, base_count: i64) -> Self {
        self.base_count = base_count;
        self
    }

    pub fn set_block_count(mut self, block_count: usize) -> Self {
        self.block_count = block_count;
        self
    }

    pub fn set_landmarks(mut self, landmarks: Vec<i32>) -> Self {
        self.landmarks = landmarks;
        self
    }

    pub fn set_crc32(mut self, crc32: u32) -> Self {
        self.crc32 = crc32;
        self
    }

    pub fn build(self) -> Header {
        Header {
            length: self.length,
            reference_sequence_id: self.reference_sequence_id,
            start_position: self.start_position,
            alignment_span: self.alignment_span,
            record_count: self.record_count,
            record_counter: self.record_counter,
            base_count: self.base_count,
            block_count: self.block_count,
            landmarks: self.landmarks,
            crc32: self.crc32,
        }
    }
}
