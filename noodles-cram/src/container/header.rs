mod builder;

pub use self::builder::Builder;

use super::ReferenceSequenceContext;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Header {
    reference_sequence_context: ReferenceSequenceContext,
    record_count: usize,
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

    pub fn reference_sequence_context(&self) -> ReferenceSequenceContext {
        self.reference_sequence_context
    }

    pub fn record_count(&self) -> usize {
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
