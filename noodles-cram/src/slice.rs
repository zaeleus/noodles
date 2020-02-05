use crate::num::Itf8;

#[derive(Debug)]
pub struct Slice {
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

impl Slice {
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
