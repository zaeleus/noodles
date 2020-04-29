use crate::num::Itf8;

#[derive(Debug)]
pub struct Header {
    reference_sequence_id: Itf8,
    alignment_start: Itf8,
    alignment_span: Itf8,
    n_records: Itf8,
    record_counter: Itf8,
    n_blocks: Itf8,
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
        n_records: Itf8,
        record_counter: Itf8,
        n_blocks: Itf8,
        block_content_ids: Vec<Itf8>,
        embedded_reference_bases_block_content_id: Itf8,
        reference_md5: [u8; 16],
        optional_tags: Vec<u8>,
    ) -> Self {
        Self {
            reference_sequence_id,
            alignment_start,
            alignment_span,
            n_records,
            record_counter,
            n_blocks,
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

    pub fn n_records(&self) -> Itf8 {
        self.n_records
    }

    pub fn n_blocks(&self) -> Itf8 {
        self.n_blocks
    }
}
