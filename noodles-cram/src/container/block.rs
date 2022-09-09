mod builder;
mod compression_method;
mod content_type;

pub use self::{
    builder::Builder, compression_method::CompressionMethod, content_type::ContentType,
};

use std::{
    io::{self, Read},
    mem,
};

use bytes::Bytes;
use xz2::read::XzDecoder;

use crate::{
    codecs::{
        aac::arith_decode, fqzcomp::fqz_decode, name_tokenizer::decode_names, rans::rans_decode,
        rans_nx16::rans_decode_nx16,
    },
    num::itf8,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Block {
    compression_method: CompressionMethod,
    content_type: ContentType,
    content_id: i32,
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

    pub fn content_id(&self) -> i32 {
        self.content_id
    }

    pub fn uncompressed_len(&self) -> usize {
        self.uncompressed_len
    }

    pub fn data(&self) -> &[u8] {
        &self.data
    }

    pub fn decompressed_data(&self) -> io::Result<Bytes> {
        use crate::codecs::{bzip2, gzip};

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
                let mut reader = XzDecoder::new(self.data());
                let mut buf = Vec::with_capacity(self.uncompressed_len);
                reader.read_to_end(&mut buf)?;
                Ok(Bytes::from(buf))
            }
            CompressionMethod::Rans4x8 => {
                let mut buf = self.data();
                rans_decode(&mut buf).map(Bytes::from)
            }
            CompressionMethod::RansNx16 => {
                let mut reader = self.data();
                rans_decode_nx16(&mut reader, self.uncompressed_len()).map(Bytes::from)
            }
            CompressionMethod::AdaptiveArithmeticCoding => {
                let mut reader = self.data();
                arith_decode(&mut reader, self.uncompressed_len()).map(Bytes::from)
            }
            CompressionMethod::Fqzcomp => {
                let mut reader = self.data();
                fqz_decode(&mut reader).map(Bytes::from)
            }
            CompressionMethod::NameTokenizer => {
                let mut reader = self.data();
                let names = decode_names(&mut reader)?;
                let data: Vec<_> = names.into_iter().flat_map(|s| s.into_bytes()).collect();
                Ok(Bytes::from(data))
            }
        }
    }

    pub fn len(&self) -> usize {
        // method
        mem::size_of::<u8>()
            // block content type ID
            + mem::size_of::<u8>()
            + itf8::size_of(self.content_id())
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
