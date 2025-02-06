mod builder;
mod compression_method;
mod content_type;

pub use self::{
    builder::Builder, compression_method::CompressionMethod, content_type::ContentType,
};

use std::{io, mem};

use bytes::Bytes;

use crate::{
    codecs::{aac, fqzcomp, name_tokenizer, rans_4x8, rans_nx16},
    num::itf8,
};

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

    pub fn decompressed_data(&self) -> io::Result<Bytes> {
        use crate::codecs::{bzip2, gzip, lzma};

        match self.compression_method {
            CompressionMethod::None => Ok(self.data.clone()),
            CompressionMethod::Gzip => {
                let mut dst = vec![0; self.uncompressed_len];
                gzip::decode(self.data(), &mut dst)?;
                Ok(Bytes::from(dst))
            }
            CompressionMethod::Bzip2 => {
                let mut dst = vec![0; self.uncompressed_len];
                bzip2::decode(self.data(), &mut dst)?;
                Ok(Bytes::from(dst))
            }
            CompressionMethod::Lzma => {
                let mut dst = vec![0; self.uncompressed_len];
                lzma::decode(self.data(), &mut dst)?;
                Ok(Bytes::from(dst))
            }
            CompressionMethod::Rans4x8 => {
                let mut buf = self.data();
                rans_4x8::decode(&mut buf).map(Bytes::from)
            }
            CompressionMethod::RansNx16 => {
                let mut reader = self.data();
                rans_nx16::decode(&mut reader, self.uncompressed_len()).map(Bytes::from)
            }
            CompressionMethod::AdaptiveArithmeticCoding => {
                let mut reader = self.data();
                aac::decode(&mut reader, self.uncompressed_len()).map(Bytes::from)
            }
            CompressionMethod::Fqzcomp => {
                let mut reader = self.data();
                fqzcomp::decode(&mut reader).map(Bytes::from)
            }
            CompressionMethod::NameTokenizer => {
                let mut reader = self.data();
                name_tokenizer::decode(&mut reader).map(Bytes::from)
            }
        }
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
