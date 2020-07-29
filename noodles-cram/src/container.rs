pub mod block;
pub mod compression_header;
mod header;
pub mod reference_sequence_id;
pub mod slice;

pub use self::{
    block::Block, compression_header::CompressionHeader, header::Header,
    reference_sequence_id::ReferenceSequenceId, slice::Slice,
};

#[derive(Debug, Default)]
pub struct Container {
    header: Header,
    blocks: Vec<Block>,
}

impl Container {
    /// Creates an EOF container.
    pub fn eof() -> Self {
        Self::new(Header::eof(), vec![Block::eof()])
    }

    pub fn new(header: Header, blocks: Vec<Block>) -> Self {
        Self { header, blocks }
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn blocks(&self) -> &[Block] {
        &self.blocks
    }

    pub fn is_eof(&self) -> bool {
        self.header.is_eof()
    }
}
