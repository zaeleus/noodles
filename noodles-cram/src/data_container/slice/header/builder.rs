use crate::{
    container::ReferenceSequenceId,
    num::{Itf8, Ltf8},
};

use super::{EmbeddedReferenceBasesBlockContentId, Header};

#[derive(Default)]
pub struct Builder {
    reference_sequence_id: ReferenceSequenceId,
    alignment_start: Itf8,
    alignment_span: Itf8,
    record_count: usize,
    record_counter: Ltf8,
    block_count: usize,
    block_content_ids: Vec<Itf8>,
    embedded_reference_bases_block_content_id: EmbeddedReferenceBasesBlockContentId,
    reference_md5: [u8; 16],
    optional_tags: Vec<u8>,
}

impl Builder {
    pub fn set_reference_sequence_id(mut self, reference_sequence_id: ReferenceSequenceId) -> Self {
        self.reference_sequence_id = reference_sequence_id;
        self
    }

    pub fn set_alignment_start(mut self, alignment_start: Itf8) -> Self {
        self.alignment_start = alignment_start;
        self
    }

    pub fn set_alignment_span(mut self, alignment_span: Itf8) -> Self {
        self.alignment_span = alignment_span;
        self
    }

    pub fn set_record_count(mut self, record_count: usize) -> Self {
        self.record_count = record_count;
        self
    }

    pub fn set_record_counter(mut self, record_counter: Ltf8) -> Self {
        self.record_counter = record_counter;
        self
    }

    pub fn set_block_count(mut self, block_count: usize) -> Self {
        self.block_count = block_count;
        self
    }

    pub fn set_block_content_ids(mut self, block_content_ids: Vec<Itf8>) -> Self {
        self.block_content_ids = block_content_ids;
        self
    }

    pub fn set_embedded_reference_bases_block_content_id(
        mut self,
        embedded_reference_bases_block_content_id: EmbeddedReferenceBasesBlockContentId,
    ) -> Self {
        self.embedded_reference_bases_block_content_id = embedded_reference_bases_block_content_id;
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
            reference_sequence_id: self.reference_sequence_id,
            alignment_start: self.alignment_start,
            alignment_span: self.alignment_span,
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
