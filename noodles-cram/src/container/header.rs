use crate::num::{Itf8, Ltf8};

use super::ReferenceSequenceId;

// § 9 End of file container (2020-06-22)
const EOF_LEN: i32 = 15;
const EOF_START_POSITION: Itf8 = 4_542_278;
const EOF_BLOCK_COUNT: Itf8 = 1;
const EOF_CRC32: u32 = 0x4f_d9_bd_05;

#[derive(Debug, Default)]
pub struct Header {
    length: i32,
    reference_sequence_id: ReferenceSequenceId,
    start_position: Itf8,
    alignment_span: Itf8,
    record_count: Itf8,
    record_counter: Ltf8,
    base_count: Ltf8,
    block_count: Itf8,
    landmarks: Vec<Itf8>,
    crc32: u32,
}

#[allow(clippy::len_without_is_empty)]
impl Header {
    /// Creates a container header used in the EOF container.
    pub fn eof() -> Self {
        Self {
            length: EOF_LEN,
            start_position: EOF_START_POSITION,
            block_count: EOF_BLOCK_COUNT,
            crc32: EOF_CRC32,
            ..Default::default()
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn new(
        length: i32,
        reference_sequence_id: ReferenceSequenceId,
        start_position: Itf8,
        alignment_span: Itf8,
        record_count: Itf8,
        record_counter: Ltf8,
        base_count: Ltf8,
        block_count: Itf8,
        landmarks: Vec<Itf8>,
        crc32: u32,
    ) -> Self {
        Self {
            length,
            reference_sequence_id,
            start_position,
            alignment_span,
            record_count,
            record_counter,
            base_count,
            block_count,
            landmarks,
            crc32,
        }
    }

    pub fn len(&self) -> i32 {
        self.length
    }

    pub fn reference_sequence_id(&self) -> ReferenceSequenceId {
        self.reference_sequence_id
    }

    pub fn start_position(&self) -> Itf8 {
        self.start_position
    }

    pub fn alignment_span(&self) -> Itf8 {
        self.alignment_span
    }

    pub fn record_count(&self) -> Itf8 {
        self.record_count
    }

    pub fn record_counter(&self) -> Ltf8 {
        self.record_counter
    }

    pub fn base_count(&self) -> Ltf8 {
        self.base_count
    }

    pub fn block_count(&self) -> Itf8 {
        self.block_count
    }

    pub fn landmarks(&self) -> &[Itf8] {
        &self.landmarks
    }

    pub fn crc32(&self) -> u32 {
        self.crc32
    }

    pub fn is_eof(&self) -> bool {
        self.length == EOF_LEN
            && self.reference_sequence_id.is_none()
            && self.start_position == EOF_START_POSITION
            && self.alignment_span == 0
            && self.record_count == 0
            && self.record_counter == 0
            && self.base_count == 0
            && self.block_count == EOF_BLOCK_COUNT
            && self.landmarks.is_empty()
            && self.crc32 == EOF_CRC32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_eof() {
        let header = Header::eof();
        assert!(header.is_eof());
    }
}
