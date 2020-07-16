pub mod block;
pub mod compression_header;
mod header;
pub mod slice;

pub use self::{block::Block, compression_header::CompressionHeader, header::Header, slice::Slice};

#[derive(Debug, Default)]
pub struct Container {
    header: Header,
    compression_header: CompressionHeader,
    slices: Vec<Slice>,
}

impl Container {
    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn header_mut(&mut self) -> &mut Header {
        &mut self.header
    }

    pub fn compression_header(&self) -> &CompressionHeader {
        &self.compression_header
    }

    pub fn compression_header_mut(&mut self) -> &mut CompressionHeader {
        &mut self.compression_header
    }

    pub fn slices(&self) -> &[Slice] {
        &self.slices
    }

    pub fn slices_mut(&mut self) -> &mut Vec<Slice> {
        &mut self.slices
    }

    pub fn add_slice(&mut self, slice: Slice) {
        self.slices.push(slice);
    }

    pub fn clear(&mut self) {
        self.header = Default::default();
        self.compression_header = Default::default();
        self.slices.clear();
    }

    pub fn is_eof(&self) -> bool {
        self.header.is_eof()
    }
}
