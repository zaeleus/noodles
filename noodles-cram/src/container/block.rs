mod builder;
mod compression_method;
mod content_type;

pub use self::{
    builder::Builder, compression_method::CompressionMethod, content_type::ContentType,
};

use std::mem;

use bytes::Bytes;

use crate::num::itf8;

pub type ContentId = i32;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Block {
    compression_method: CompressionMethod,
    content_type: ContentType,
    content_id: ContentId,
    uncompressed_len: usize,
    data: Bytes,
}

#[allow(clippy::len_without_is_empty)]
impl Block {
    pub fn builder() -> Builder {
        Builder::default()
    }

    pub fn compression_method(&self) -> CompressionMethod {
        self.compression_method
    }

    pub fn content_type(&self) -> ContentType {
        self.content_type
    }

    pub fn content_id(&self) -> ContentId {
        self.content_id
    }

    pub fn uncompressed_len(&self) -> usize {
        self.uncompressed_len
    }

    pub fn data(&self) -> &[u8] {
        &self.data
    }

    pub fn len(&self) -> usize {
        // method
        mem::size_of::<u8>()
            // block content type ID
            + mem::size_of::<u8>()
            + itf8::size_of(self.content_id)
            + itf8::size_of(self.data.len() as i32)
            + itf8::size_of(self.uncompressed_len() as i32)
            + self.data.len()
            // crc32
            + mem::size_of::<u32>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_len() {
        let data = Bytes::from_static(b"noodles");

        let block = Block::builder()
            .set_content_type(ContentType::ExternalData)
            .set_uncompressed_len(data.len())
            .set_data(data)
            .build();

        assert_eq!(block.len(), 16);
    }
}
