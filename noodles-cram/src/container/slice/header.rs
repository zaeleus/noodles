use crate::container::{block, ReferenceSequenceContext};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Header {
    pub(crate) reference_sequence_context: ReferenceSequenceContext,
    pub(crate) record_count: usize,
    pub(crate) record_counter: u64,
    pub(crate) block_count: usize,
    pub(crate) block_content_ids: Vec<block::ContentId>,
    pub(crate) embedded_reference_bases_block_content_id: Option<block::ContentId>,
    pub(crate) reference_md5: [u8; 16],
    pub(crate) optional_tags: Vec<u8>,
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

    pub fn block_count(&self) -> usize {
        self.block_count
    }

    pub fn block_content_ids(&self) -> &[block::ContentId] {
        &self.block_content_ids
    }

    pub fn embedded_reference_bases_block_content_id(&self) -> Option<block::ContentId> {
        self.embedded_reference_bases_block_content_id
    }

    pub fn reference_md5(&self) -> &[u8] {
        &self.reference_md5
    }

    pub fn optional_tags(&self) -> &[u8] {
        &self.optional_tags
    }
}
