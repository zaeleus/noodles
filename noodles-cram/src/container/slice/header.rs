mod builder;
mod embedded_reference_bases_block_content_id;

pub use {
    builder::Builder,
    embedded_reference_bases_block_content_id::EmbeddedReferenceBasesBlockContentId,
};

use crate::{
    container::ReferenceSequenceId,
    num::{Itf8, Ltf8},
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Header {
    reference_sequence_id: ReferenceSequenceId,
    alignment_start: Itf8,
    alignment_span: Itf8,
    record_count: Itf8,
    record_counter: Ltf8,
    block_count: Itf8,
    block_content_ids: Vec<Itf8>,
    embedded_reference_bases_block_content_id: EmbeddedReferenceBasesBlockContentId,
    reference_md5: [u8; 16],
    optional_tags: Vec<u8>,
}

impl Header {
    pub fn builder() -> Builder {
        Builder::default()
    }

    pub fn reference_sequence_id(&self) -> ReferenceSequenceId {
        self.reference_sequence_id
    }

    pub fn alignment_start(&self) -> Itf8 {
        self.alignment_start
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

    pub fn block_count(&self) -> Itf8 {
        self.block_count
    }

    pub fn block_content_ids(&self) -> &[Itf8] {
        &self.block_content_ids
    }

    pub fn embedded_reference_bases_block_content_id(
        &self,
    ) -> EmbeddedReferenceBasesBlockContentId {
        self.embedded_reference_bases_block_content_id
    }

    pub fn reference_md5(&self) -> &[u8] {
        &self.reference_md5
    }

    pub fn optional_tags(&self) -> &[u8] {
        &self.optional_tags
    }
}
