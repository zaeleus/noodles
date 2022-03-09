mod builder;

pub use self::builder::Builder;

use noodles_sam as sam;

use super::ReferenceSequenceId;

// § 9 End of file container (2020-06-22)
const EOF_LEN: i32 = 15;
const EOF_START_POSITION: i32 = 4_542_278;
const EOF_BLOCK_COUNT: usize = 1;
const EOF_CRC32: u32 = 0x4f_d9_bd_05;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Header {
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

#[allow(clippy::len_without_is_empty)]
impl Header {
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Creates a container header used in the EOF container.
    pub fn eof() -> Self {
        Self {
            length: EOF_LEN,
            start_position: Some(sam::record::Position::try_from(EOF_START_POSITION).unwrap()),
            block_count: EOF_BLOCK_COUNT,
            crc32: EOF_CRC32,
            ..Default::default()
        }
    }

    pub fn len(&self) -> i32 {
        self.length
    }

    pub fn reference_sequence_id(&self) -> ReferenceSequenceId {
        self.reference_sequence_id
    }

    pub fn start_position(&self) -> Option<sam::record::Position> {
        self.start_position
    }

    pub fn alignment_span(&self) -> usize {
        self.alignment_span
    }

    pub fn record_count(&self) -> i32 {
        self.record_count
    }

    pub fn record_counter(&self) -> i64 {
        self.record_counter
    }

    pub fn base_count(&self) -> i64 {
        self.base_count
    }

    pub fn block_count(&self) -> usize {
        self.block_count
    }

    pub fn landmarks(&self) -> &[i32] {
        &self.landmarks
    }

    pub fn is_eof(&self) -> bool {
        self.length == EOF_LEN
            && self.reference_sequence_id.is_none()
            && self
                .start_position
                .map(|position| i32::from(position) == EOF_START_POSITION)
                .unwrap_or(false)
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
