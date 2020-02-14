use crate::{num::Itf8, reader::record, Block, CompressionHeader};

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
    pub fn n_blocks(&self) -> Itf8 {
        self.n_blocks
    }
}

#[derive(Debug)]
pub struct Slice {
    header: Header,
    core_block: Block,
    external_blocks: Vec<Block>,
}

impl Slice {
    pub fn new(header: Header) -> Self {
        Self {
            header,
            core_block: Block::default(),
            external_blocks: Vec::new(),
        }
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn core_block_mut(&mut self) -> &mut Block {
        &mut self.core_block
    }

    pub fn add_external_block(&mut self, block: Block) {
        self.external_blocks.push(block);
    }

    pub fn records<'a>(&'a self, compression_header: &'a CompressionHeader) -> record::Reader<'a> {
        record::Reader::new(compression_header, &self.external_blocks)
    }
}
