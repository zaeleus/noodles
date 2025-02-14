mod compression_method;
mod content_type;

use std::{
    io::{self, Write},
    mem,
};

use byteorder::{LittleEndian, WriteBytesExt};
use flate2::CrcWriter;

use self::{compression_method::write_compression_method, content_type::write_content_type};
use crate::{
    codecs::Encoder,
    container::block::{CompressionMethod, ContentId, ContentType},
    io::writer::num::write_itf8,
};

pub struct Block {
    pub(crate) compression_method: CompressionMethod,
    pub(crate) content_type: ContentType,
    pub(crate) content_id: ContentId,
    pub(crate) uncompressed_size: usize,
    pub(crate) src: Vec<u8>,
}

impl Block {
    pub fn encode(
        content_type: ContentType,
        content_id: ContentId,
        encoder: Option<&Encoder>,
        src: &[u8],
    ) -> io::Result<Self> {
        use crate::codecs::{aac, bzip2, gzip, lzma, name_tokenizer, rans_4x8, rans_nx16};

        let (compression_method, buf) = match encoder {
            None => (CompressionMethod::None, src.to_vec()),
            Some(Encoder::Gzip(compression_level)) => (
                CompressionMethod::Gzip,
                gzip::encode(*compression_level, src)?,
            ),
            Some(Encoder::Bzip2(compression_level)) => (
                CompressionMethod::Bzip2,
                bzip2::encode(*compression_level, src)?,
            ),
            Some(Encoder::Lzma(compression_level)) => (
                CompressionMethod::Lzma,
                lzma::encode(*compression_level, src)?,
            ),
            Some(Encoder::Rans4x8(order)) => {
                (CompressionMethod::Rans4x8, rans_4x8::encode(*order, src)?)
            }
            Some(Encoder::RansNx16(flags)) => {
                (CompressionMethod::RansNx16, rans_nx16::encode(*flags, src)?)
            }
            Some(Encoder::AdaptiveArithmeticCoding(flags)) => (
                CompressionMethod::AdaptiveArithmeticCoding,
                aac::encode(*flags, src)?,
            ),
            Some(Encoder::NameTokenizer) => (
                CompressionMethod::NameTokenizer,
                name_tokenizer::encode(src)?,
            ),
            Some(Encoder::Fqzcomp) => unimplemented!(),
        };

        Ok(Self {
            compression_method,
            content_type,
            content_id,
            uncompressed_size: src.len(),
            src: buf,
        })
    }

    pub fn size(&self) -> io::Result<usize> {
        let compressed_size = i32::try_from(self.src.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let uncompressed_size = i32::try_from(self.uncompressed_size)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        Ok(
            mem::size_of::<u8>() // method
                + mem::size_of::<u8>() // block content type ID
                + itf8_size_of(self.content_id)
                + itf8_size_of(compressed_size)
                + itf8_size_of(uncompressed_size)
                + self.src.len()
                + mem::size_of::<u32>(), // CRC32
        )
    }
}

pub fn write_block<W>(writer: &mut W, block: &Block) -> io::Result<()>
where
    W: Write,
{
    let mut crc_writer = CrcWriter::new(writer);
    write_block_inner(&mut crc_writer, block)
}

fn write_block_inner<W>(writer: &mut CrcWriter<W>, block: &Block) -> io::Result<()>
where
    W: Write,
{
    write_compression_method(writer, block.compression_method)?;

    write_content_type(writer, block.content_type)?;
    write_itf8(writer, block.content_id)?;

    write_size(writer, block.src.len())?; // compressed size
    write_size(writer, block.uncompressed_size)?;

    writer.write_all(&block.src)?;

    let crc32 = writer.crc().sum();
    writer.get_mut().write_u32::<LittleEndian>(crc32)?;

    Ok(())
}

fn write_size<W>(writer: &mut W, size: usize) -> io::Result<()>
where
    W: Write,
{
    let n = i32::try_from(size).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(writer, n)
}

pub fn itf8_size_of(n: i32) -> usize {
    if n >> (8 - 1) == 0 {
        1
    } else if n >> (16 - 2) == 0 {
        2
    } else if n >> (24 - 3) == 0 {
        3
    } else if n >> (32 - 4) == 0 {
        4
    } else {
        5
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::container::block::{CompressionMethod, ContentType};

    #[test]
    fn test_write_block() -> io::Result<()> {
        let block = Block {
            compression_method: CompressionMethod::None,
            content_type: ContentType::ExternalData,
            content_id: 1,
            uncompressed_size: 4,
            src: b"ndls".to_vec(),
        };

        let mut buf = Vec::new();
        write_block(&mut buf, &block)?;

        let expected = [
            0x00, // compression method = none
            0x04, // content type = external data
            0x01, // content ID = 1
            0x04, // uncompressed size = 4
            0x04, // compressed size = 4
            b'n', b'd', b'l', b's', // data = b"ndls"
            0xd7, 0x12, 0x46, 0x3e, // CRC32 = 3e4612d7
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_itf8_size_of() {
        assert_eq!(itf8_size_of(0), 1);
        assert_eq!(itf8_size_of(1877), 2);
        assert_eq!(itf8_size_of(480665), 3);
        assert_eq!(itf8_size_of(123050342), 4);
        assert_eq!(itf8_size_of(1968805474), 5);
        assert_eq!(itf8_size_of(-1), 5);
    }
}
