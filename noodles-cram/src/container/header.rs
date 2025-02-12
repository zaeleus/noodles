use super::ReferenceSequenceContext;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Header {
    pub reference_sequence_context: ReferenceSequenceContext,
    pub record_count: usize,
    pub record_counter: u64,
    pub base_count: u64,
    pub block_count: usize,
    pub landmarks: Vec<usize>,
}

impl Header {
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
