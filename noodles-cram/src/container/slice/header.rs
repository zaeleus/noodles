use crate::num::{Itf8, Ltf8};

#[derive(Debug)]
pub struct Header {
    reference_sequence_id: Itf8,
    alignment_start: Itf8,
    alignment_span: Itf8,
    record_count: Itf8,
    record_counter: Ltf8,
    block_count: Itf8,
    block_content_ids: Vec<Itf8>,
    embedded_reference_bases_block_content_id: Itf8,
    reference_md5: [u8; 16],
    optional_tags: Vec<u8>,
}

impl Header {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        reference_sequence_id: Itf8,
        alignment_start: Itf8,
        alignment_span: Itf8,
        record_count: Itf8,
        record_counter: Ltf8,
        block_count: Itf8,
        block_content_ids: Vec<Itf8>,
        embedded_reference_bases_block_content_id: Itf8,
        reference_md5: [u8; 16],
        optional_tags: Vec<u8>,
    ) -> Self {
        Self {
            reference_sequence_id,
            alignment_start,
            alignment_span,
            record_count,
            record_counter,
            block_count,
            block_content_ids,
            embedded_reference_bases_block_content_id,
            reference_md5,
            optional_tags,
        }
    }
}

impl Header {
    pub fn reference_sequence_id(&self) -> Itf8 {
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

    pub fn embedded_reference_bases_block_content_id(&self) -> Itf8 {
        self.embedded_reference_bases_block_content_id
    }

    pub fn reference_md5(&self) -> &[u8] {
        &self.reference_md5
    }

    pub fn optional_tags(&self) -> &[u8] {
        &self.optional_tags
    }
}
