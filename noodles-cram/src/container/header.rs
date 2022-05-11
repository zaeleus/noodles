mod builder;

pub use self::builder::Builder;

use noodles_core::Position;

use super::ReferenceSequenceId;

// § 9 End of file container (2020-06-22)
const EOF_LEN: usize = 15;
const EOF_START_POSITION: usize = 4_542_278;
const EOF_BLOCK_COUNT: usize = 1;
const EOF_CRC32: u32 = 0x4f_d9_bd_05;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Header {
    length: usize,
    reference_sequence_id: ReferenceSequenceId,
    start_position: Option<Position>,
    alignment_span: usize,
    record_count: i32,
    record_counter: i64,
    base_count: u64,
    block_count: usize,
    landmarks: Vec<usize>,
    crc32: u32,
}

#[allow(clippy::len_without_is_empty)]
impl Header {
    pub fn builder() -> Builder {
        Builder::default()
    }

    pub fn len(&self) -> usize {
        self.length
    }

    pub fn reference_sequence_id(&self) -> ReferenceSequenceId {
        self.reference_sequence_id
    }

    pub fn start_position(&self) -> Option<Position> {
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

    pub fn base_count(&self) -> u64 {
        self.base_count
    }

    pub fn block_count(&self) -> usize {
        self.block_count
    }

    pub fn landmarks(&self) -> &[usize] {
        &self.landmarks
    }

    pub fn is_eof(&self) -> bool {
        self.length == EOF_LEN
            && self.reference_sequence_id.is_none()
            && self
                .start_position
                .map(|position| usize::from(position) == EOF_START_POSITION)
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
