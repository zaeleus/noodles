pub(crate) mod compression_method;
pub(crate) mod content_type;

use std::io;

use flate2::Crc;

use self::{compression_method::read_compression_method, content_type::read_content_type};
use crate::{
    container::block::{CompressionMethod, ContentId, ContentType},
    io::reader::num::{read_itf8, read_itf8_as, read_u32_le},
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Block<'c> {
    pub compression_method: CompressionMethod,
    pub content_type: ContentType,
    pub content_id: ContentId,
    pub uncompressed_size: usize,
    pub src: &'c [u8],
}

impl Block<'_> {
    pub fn decode(&self) -> io::Result<Vec<u8>> {
        use crate::codecs::{aac, bzip2, fqzcomp, gzip, lzma, name_tokenizer, rans_4x8, rans_nx16};

        match self.compression_method {
            CompressionMethod::None => Ok(self.src.to_vec()),
            CompressionMethod::Gzip => {
                let mut dst = vec![0; self.uncompressed_size];
                gzip::decode(self.src, &mut dst)?;
                Ok(dst)
            }
            CompressionMethod::Bzip2 => {
                let mut dst = vec![0; self.uncompressed_size];
                bzip2::decode(self.src, &mut dst)?;
                Ok(dst)
            }
            CompressionMethod::Lzma => {
                let mut dst = vec![0; self.uncompressed_size];
                lzma::decode(self.src, &mut dst)?;
                Ok(dst)
            }
            CompressionMethod::Rans4x8 => rans_4x8::decode(&mut &self.src[..]),
            CompressionMethod::RansNx16 => rans_nx16::decode(self.src, self.uncompressed_size),
            CompressionMethod::AdaptiveArithmeticCoding => {
                aac::decode(&mut &self.src[..], self.uncompressed_size)
            }
            CompressionMethod::Fqzcomp => fqzcomp::decode(&mut &self.src[..]),
            CompressionMethod::NameTokenizer => name_tokenizer::decode(self.src),
        }
    }
}

fn read_block<'c>(src: &mut &'c [u8]) -> io::Result<Block<'c>> {
    let original_src = *src;

    let mut compression_method = read_compression_method(src)?;

    let content_type = read_content_type(src)?;
    let content_id = read_itf8(src)?;

    let compressed_size = read_itf8_as(src)?;
    let uncompressed_size = read_itf8_as(src)?;

    let (data, rest) = src
        .split_at_checked(compressed_size)
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    let end = original_src.len() - src.len();
    let actual_crc32 = crc32(&original_src[..end]);

    let expected_crc32 = read_u32_le(src)?;

    if actual_crc32 != expected_crc32 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "container block checksum mismatch: expected {expected_crc32:08x}, got {actual_crc32:08x}"
            ),
        ));
    }

    // § 8 "Block structure" (2024-09-04): "Blocks with a raw (uncompressed) size of zero are
    // treated as empty, irrespective of their `method` byte."
    if uncompressed_size == 0 {
        compression_method = CompressionMethod::None;
    }

    Ok(Block {
        compression_method,
        content_type,
        content_id,
        uncompressed_size,
        src: data,
    })
}

pub fn read_block_as<'c>(src: &mut &'c [u8], content_type: ContentType) -> io::Result<Block<'c>> {
    let block = read_block(src)?;
    validate_content_type(block.content_type, content_type)?;
    Ok(block)
}

fn crc32(buf: &[u8]) -> u32 {
    let mut crc = Crc::new();
    crc.update(buf);
    crc.sum()
}

fn validate_content_type(actual: ContentType, expected: ContentType) -> io::Result<()> {
    if actual == expected {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid block content type: expected {expected:?}, got {actual:?}"),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_block() -> io::Result<()> {
        let src = [
            0x00, // compression method = none (0)
            0x04, // content type = external data (4)
            0x01, // block content ID = 1
            0x04, // size in bytes = 4 bytes
            0x04, // raw size in bytes = 4 bytes
            0x6e, 0x64, 0x6c, 0x73, // data = b"ndls",
            0xd7, 0x12, 0x46, 0x3e, // CRC32 = 3e4612d7
        ];

        let actual = read_block(&mut &src[..])?;

        let expected = Block {
            compression_method: CompressionMethod::None,
            content_type: ContentType::ExternalData,
            content_id: ContentId::from(1),
            uncompressed_size: 4,
            src: b"ndls",
        };

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_block_with_empty_block() -> io::Result<()> {
        let src = [
            0x04, // compression method = rANS 4x8 (4)
            0x04, // content type = external data (4)
            0x01, // block content ID = 1
            0x00, // size in bytes = 0 bytes
            0x00, // raw size in bytes = 0 bytes
            // data = b"",
            0xbd, 0xac, 0x02, 0xbd, // CRC32 = bd02acbd
        ];

        let actual = read_block(&mut &src[..])?;

        let expected = Block {
            compression_method: CompressionMethod::None,
            content_type: ContentType::ExternalData,
            content_id: ContentId::from(1),
            uncompressed_size: 0,
            src: &[],
        };

        assert_eq!(actual, expected);

        Ok(())
    }
}
