use crate::num::Itf8;

use super::{Block, CompressionMethod, ContentType};

#[derive(Debug, Default)]
pub struct Builder {
    compression_method: CompressionMethod,
    content_type: Option<ContentType>,
    content_id: Itf8,
    uncompressed_len: Itf8,
    data: Vec<u8>,
    crc32: u32,
}

impl Builder {
    pub fn set_compression_method(mut self, compression_method: CompressionMethod) -> Self {
        self.compression_method = compression_method;
        self
    }

    pub fn set_content_type(mut self, content_type: ContentType) -> Self {
        self.content_type = Some(content_type);
        self
    }

    pub fn set_content_id(mut self, content_id: Itf8) -> Self {
        self.content_id = content_id;
        self
    }

    pub fn set_uncompressed_len(mut self, uncompressed_len: Itf8) -> Self {
        self.uncompressed_len = uncompressed_len;
        self
    }

    pub fn set_data(mut self, data: Vec<u8>) -> Self {
        self.data = data;
        self
    }

    pub fn set_crc32(mut self, crc32: u32) -> Self {
        self.crc32 = crc32;
        self
    }

    pub fn build(self) -> Block {
        Block {
            compression_method: self.compression_method,
            content_type: self.content_type.expect("missing block content type"),
            content_id: self.content_id,
            uncompressed_len: self.uncompressed_len,
            data: self.data,
            crc32: self.crc32,
        }
    }
}
