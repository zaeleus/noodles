use super::Header;
use crate::container::{block, ReferenceSequenceContext};

#[derive(Default)]
pub struct Builder {
    reference_sequence_context: ReferenceSequenceContext,
    record_count: usize,
    record_counter: u64,
    block_count: usize,
    block_content_ids: Vec<block::ContentId>,
    embedded_reference_bases_block_content_id: Option<block::ContentId>,
    reference_md5: [u8; 16],
    optional_tags: Vec<u8>,
}

impl Builder {
    pub fn set_reference_sequence_context(
        mut self,
        reference_sequence_context: ReferenceSequenceContext,
    ) -> Self {
        self.reference_sequence_context = reference_sequence_context;
        self
    }

    pub fn set_record_count(mut self, record_count: usize) -> Self {
        self.record_count = record_count;
        self
    }

    pub fn set_record_counter(mut self, record_counter: u64) -> Self {
        self.record_counter = record_counter;
        self
    }

    pub fn set_block_count(mut self, block_count: usize) -> Self {
        self.block_count = block_count;
        self
    }

    pub fn set_block_content_ids(mut self, block_content_ids: Vec<block::ContentId>) -> Self {
        self.block_content_ids = block_content_ids;
        self
    }

    pub fn set_embedded_reference_bases_block_content_id(mut self, id: block::ContentId) -> Self {
        self.embedded_reference_bases_block_content_id = Some(id);
        self
    }

    pub fn set_reference_md5(mut self, reference_md5: [u8; 16]) -> Self {
        self.reference_md5 = reference_md5;
        self
    }

    pub fn set_optional_tags(mut self, optional_tags: Vec<u8>) -> Self {
        self.optional_tags = optional_tags;
        self
    }

    pub fn build(self) -> Header {
        Header {
            reference_sequence_context: self.reference_sequence_context,
            record_count: self.record_count,
            record_counter: self.record_counter,
            block_count: self.block_count,
            block_content_ids: self.block_content_ids,
            embedded_reference_bases_block_content_id: self
                .embedded_reference_bases_block_content_id,
            reference_md5: self.reference_md5,
            optional_tags: self.optional_tags,
        }
    }
}
