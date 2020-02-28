use std::{collections::HashMap, io::Cursor};

use crate::{num::Itf8, reader::record, BitReader, Block, CompressionHeader};

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
    pub fn n_records(&self) -> Itf8 {
        self.n_records
    }

    pub fn n_blocks(&self) -> Itf8 {
        self.n_blocks
    }
}

#[derive(Debug)]
pub struct Slice {
    header: Header,
    core_data_block: Block,
    external_blocks: Vec<Block>,
}

impl Slice {
    pub fn new(header: Header) -> Self {
        Self {
            header,
            core_data_block: Block::default(),
            external_blocks: Vec::new(),
        }
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn core_data_block_mut(&mut self) -> &mut Block {
        &mut self.core_data_block
    }

    pub fn add_external_block(&mut self, block: Block) {
        self.external_blocks.push(block);
    }

    pub fn records<'a>(
        &'a self,
        compression_header: &'a CompressionHeader,
    ) -> record::Reader<'a, Cursor<Vec<u8>>, Cursor<Vec<u8>>> {
        let core_data = BitReader::new(Cursor::new(self.core_data_block.decompressed_data()));

        let external_data: HashMap<_, _> = self
            .external_blocks
            .iter()
            .map(|block| (block.content_id(), Cursor::new(block.decompressed_data())))
            .collect();

        record::Reader::new(
            compression_header,
            core_data,
            external_data,
            self.header.reference_sequence_id,
            self.header.alignment_start,
        )
    }
}
