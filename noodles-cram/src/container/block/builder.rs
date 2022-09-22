use std::io;

use bytes::Bytes;

use super::{Block, CompressionMethod, ContentId, ContentType};
use crate::codecs::Encoder;

#[derive(Debug, Default)]
pub struct Builder {
    compression_method: CompressionMethod,
    content_type: Option<ContentType>,
    content_id: ContentId,
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

    pub fn set_content_id(mut self, content_id: ContentId) -> Self {
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
    pub fn compress_and_set_data(mut self, data: Vec<u8>, encoder: Encoder) -> io::Result<Self> {
        use crate::codecs::{aac, bzip2, gzip, lzma, name_tokenizer, rans_4x8, rans_nx16};

        self.uncompressed_len = data.len();

        let (compression_method, data) = match encoder {
            Encoder::Gzip(compression_level) => (
                CompressionMethod::Gzip,
                gzip::encode(compression_level, &data)?,
            ),
            Encoder::Bzip2(compression_level) => (
                CompressionMethod::Bzip2,
                bzip2::encode(compression_level, &data)?,
            ),
            Encoder::Lzma(compression_level) => (
                CompressionMethod::Lzma,
                lzma::encode(compression_level, &data)?,
            ),
            Encoder::Rans4x8(order) => {
                (CompressionMethod::Rans4x8, rans_4x8::encode(order, &data)?)
            }
            Encoder::RansNx16(flags) => (
                CompressionMethod::RansNx16,
                rans_nx16::encode(flags, &data)?,
            ),
            Encoder::AdaptiveArithmeticCoding(flags) => (
                CompressionMethod::AdaptiveArithmeticCoding,
                aac::encode(flags, &data)?,
            ),
            Encoder::NameTokenizer => (
                CompressionMethod::NameTokenizer,
                name_tokenizer::encode(&data)?,
            ),
        };

        self.compression_method = compression_method;
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
