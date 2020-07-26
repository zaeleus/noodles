mod header;

pub use self::header::Header;

use std::{
    borrow::Cow,
    collections::HashMap,
    convert::TryFrom,
    io::{self, Cursor},
};

use crate::{reader, BitReader};

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

    #[allow(clippy::type_complexity)]
    pub fn records<'a>(
        &'a self,
        compression_header: &'a CompressionHeader,
    ) -> reader::record::Reader<'a, Cursor<Cow<[u8]>>, Cursor<Cow<[u8]>>> {
        let core_data = BitReader::new(Cursor::new(self.core_data_block.decompressed_data()));

        let external_data: HashMap<_, _> = self
            .external_blocks
            .iter()
            .map(|block| (block.content_id(), Cursor::new(block.decompressed_data())))
            .collect();

        reader::record::Reader::new(
            compression_header,
            core_data,
            external_data,
            self.header.reference_sequence_id(),
            self.header.alignment_start(),
        )
    }
}

impl TryFrom<&[Block]> for Slice {
    type Error = io::Error;

    fn try_from(blocks: &[Block]) -> Result<Self, Self::Error> {
        let data = blocks[0].decompressed_data();
        let mut reader = &data[..];
        let header = reader::slice::read_header(&mut reader)?;

        let core_data_block = blocks[1].clone();

        let external_blocks_len = header.block_count() as usize;
        let mut external_blocks = Vec::with_capacity(external_blocks_len);

        for i in 0..external_blocks_len {
            external_blocks.push(blocks[i + 1].clone());
        }

        Ok(Slice {
            header,
            core_data_block,
            external_blocks,
        })
    }
}
