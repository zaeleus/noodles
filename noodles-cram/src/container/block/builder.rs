use std::io;

use super::{Block, CompressionMethod, ContentType};

use bytes::Bytes;

#[derive(Debug, Default)]
pub struct Builder {
    compression_method: CompressionMethod,
    content_type: Option<ContentType>,
    content_id: i32,
    uncompressed_len: usize,
    data: Bytes,
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

    pub fn set_content_id(mut self, content_id: i32) -> Self {
        self.content_id = content_id;
        self
    }

    pub fn set_uncompressed_len(mut self, uncompressed_len: usize) -> Self {
        self.uncompressed_len = uncompressed_len;
        self
    }

    pub fn set_data(mut self, data: Bytes) -> Self {
        self.data = data;
        self
    }

    /// Compresses the given data using the given compression method.
    ///
    /// This sets the compression method, the uncompressed size to the length of the given data,
    /// and the data to the compressed output of the given data.
    pub fn compress_and_set_data(
        mut self,
        data: Vec<u8>,
        compression_method: CompressionMethod,
    ) -> io::Result<Self> {
        use crate::codecs::{bzip2, gzip, lzma};

        self.compression_method = compression_method;
        self.uncompressed_len = data.len();

        let data = match compression_method {
            CompressionMethod::None => data,
            CompressionMethod::Gzip => gzip::encode(&data)?,
            CompressionMethod::Bzip2 => bzip2::encode(&data)?,
            CompressionMethod::Lzma => lzma::encode(&data)?,
            _ => unimplemented!(
                "compress_and_set_data: unhandled compression method: {:?}",
                compression_method
            ),
        };

        self.data = Bytes::from(data);

        Ok(self)
    }

    pub fn build(self) -> Block {
        Block {
            compression_method: self.compression_method,
            content_type: self.content_type.expect("missing block content type"),
            content_id: self.content_id,
            uncompressed_len: self.uncompressed_len,
            data: self.data,
        }
    }
}
