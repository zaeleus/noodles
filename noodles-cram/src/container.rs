pub mod block;
pub mod compression_header;
mod header;
pub mod slice;

pub use self::{block::Block, compression_header::CompressionHeader, header::Header, slice::Slice};

#[derive(Debug, Default)]
pub struct Container {
    header: Header,
    blocks: Vec<Block>,
}

impl Container {
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
