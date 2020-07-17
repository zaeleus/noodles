mod header;

pub use self::header::Header;

use std::{borrow::Cow, collections::HashMap, io::Cursor};

use crate::{reader::record, BitReader};

use super::{Block, CompressionHeader};

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
    ) -> record::Reader<'a, Cursor<Cow<[u8]>>, Cursor<Cow<[u8]>>> {
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
            self.header.reference_sequence_id(),
            self.header.alignment_start(),
        )
    }
}
