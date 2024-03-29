mod builder;

pub use self::builder::Builder;

use super::ReferenceSequenceContext;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Header {
    length: usize,
    reference_sequence_context: ReferenceSequenceContext,
    record_count: i32,
    record_counter: u64,
    base_count: u64,
    block_count: usize,
    landmarks: Vec<usize>,
}

#[allow(clippy::len_without_is_empty)]
impl Header {
    pub fn builder() -> Builder {
        Builder::default()
    }

    pub fn len(&self) -> usize {
        self.length
    }

    pub fn reference_sequence_context(&self) -> ReferenceSequenceContext {
        self.reference_sequence_context
    }

    pub fn record_count(&self) -> i32 {
        self.record_count
    }

    pub fn record_counter(&self) -> u64 {
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
}
