mod compression_method;
mod content_type;

use std::{
    io::{self, Write},
    mem,
};

use flate2::CrcWriter;

use self::{compression_method::write_compression_method, content_type::write_content_type};
use crate::{
    codecs::Encoder,
    container::block::{CompressionMethod, ContentId, ContentType},
    file_definition::Version,
    io::writer::num::{int_size_of, write_int, write_u32_le},
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
            Some(Encoder::Fqzcomp) => {
                return Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "fqzcomp requires record lengths; use the slice-level encoder instead",
                ));
            }
        };

        Ok(Self {
            compression_method,
            content_type,
            content_id,
            uncompressed_size: src.len(),
            src: buf,
        })
    }

    pub fn size(&self, version: Version) -> io::Result<usize> {
        let compressed_size = i32::try_from(self.src.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let uncompressed_size = i32::try_from(self.uncompressed_size)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let mut size = mem::size_of::<u8>() // method
            + mem::size_of::<u8>() // block content type ID
            + int_size_of(version, self.content_id)
            + int_size_of(version, compressed_size)
            + int_size_of(version, uncompressed_size)
            + self.src.len();

        if version.has_crc32() {
            size += mem::size_of::<u32>();
        }

        Ok(size)
    }
}

pub fn write_block<W>(writer: &mut W, block: &Block, version: Version) -> io::Result<()>
where
    W: Write,
{
    if version.has_crc32() {
        let mut crc_writer = CrcWriter::new(writer);
        write_block_body(&mut crc_writer, block, version)?;
        let crc32 = crc_writer.crc().sum();
        write_u32_le(crc_writer.get_mut(), crc32)?;
        Ok(())
    } else {
        write_block_body(writer, block, version)
    }
}

fn write_block_body<W>(writer: &mut W, block: &Block, version: Version) -> io::Result<()>
where
    W: Write,
{
    write_compression_method(writer, block.compression_method)?;

    write_content_type(writer, block.content_type)?;
    write_int(writer, version, block.content_id)?;

    write_size(writer, block.src.len(), version)?; // compressed size
    write_size(writer, block.uncompressed_size, version)?;

    writer.write_all(&block.src)?;

    Ok(())
}

fn write_size<W>(writer: &mut W, size: usize, version: Version) -> io::Result<()>
where
    W: Write,
{
    let n = i32::try_from(size).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_int(writer, version, n)
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
        write_block(&mut buf, &block, Version::default())?;

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
    fn test_int_size_of() {
        // ITF8 sizes (default version is 3.0)
        let v3 = Version::default();
        assert_eq!(int_size_of(v3, 0), 1);
        assert_eq!(int_size_of(v3, 1877), 2);
        assert_eq!(int_size_of(v3, 480665), 3);
        assert_eq!(int_size_of(v3, 123050342), 4);
        assert_eq!(int_size_of(v3, 1968805474), 5);
        assert_eq!(int_size_of(v3, -1), 5);

        // VLQ (uint7) sizes for CRAM 4.0
        let v4 = Version::new(4, 0);
        assert_eq!(int_size_of(v4, 0), 1);
        assert_eq!(int_size_of(v4, 127), 1);
        assert_eq!(int_size_of(v4, 128), 2);
        assert_eq!(int_size_of(v4, 16383), 2);
        assert_eq!(int_size_of(v4, 16384), 3);
    }
}
