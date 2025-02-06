use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use flate2::CrcWriter;

use crate::{
    container::{
        block::{CompressionMethod, ContentType},
        Block,
    },
    io::writer::num::write_itf8,
};

pub fn write_block<W>(writer: &mut W, block: &Block) -> io::Result<()>
where
    W: Write,
{
    let mut crc_writer = CrcWriter::new(writer);

    write_compression_method(&mut crc_writer, block.compression_method())?;
    write_content_type(&mut crc_writer, block.content_type())?;
    write_itf8(&mut crc_writer, block.content_id())?;

    let size_in_bytes = i32::try_from(block.data().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut crc_writer, size_in_bytes)?;

    let uncompressed_data_len = i32::try_from(block.uncompressed_len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_itf8(&mut crc_writer, uncompressed_data_len)?;

    crc_writer.write_all(block.data())?;

    let crc32 = crc_writer.crc().sum();
    let writer = crc_writer.into_inner();
    writer.write_u32::<LittleEndian>(crc32)?;

    Ok(())
}

fn write_compression_method<W>(
    writer: &mut W,
    compression_method: CompressionMethod,
) -> io::Result<()>
where
    W: Write,
{
    let n = match compression_method {
        CompressionMethod::None => 0,
        CompressionMethod::Gzip => 1,
        CompressionMethod::Bzip2 => 2,
        CompressionMethod::Lzma => 3,
        CompressionMethod::Rans4x8 => 4,
        CompressionMethod::RansNx16 => 5,
        CompressionMethod::AdaptiveArithmeticCoding => 6,
        CompressionMethod::Fqzcomp => 7,
        CompressionMethod::NameTokenizer => 8,
    };

    writer.write_u8(n)
}

fn write_content_type<W>(writer: &mut W, content_type: ContentType) -> io::Result<()>
where
    W: Write,
{
    let n = match content_type {
        ContentType::FileHeader => 0,
        ContentType::CompressionHeader => 1,
        ContentType::SliceHeader => 2,
        ContentType::Reserved => 3,
        ContentType::ExternalData => 4,
        ContentType::CoreData => 5,
    };

    writer.write_u8(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_compression_method() -> io::Result<()> {
        fn t(
            buf: &mut Vec<u8>,
            compression_method: CompressionMethod,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_compression_method(buf, compression_method)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, CompressionMethod::None, &[0x00])?;
        t(&mut buf, CompressionMethod::Gzip, &[0x01])?;
        t(&mut buf, CompressionMethod::Bzip2, &[0x02])?;
        t(&mut buf, CompressionMethod::Lzma, &[0x03])?;
        t(&mut buf, CompressionMethod::Rans4x8, &[0x04])?;
        t(&mut buf, CompressionMethod::RansNx16, &[0x05])?;
        t(
            &mut buf,
            CompressionMethod::AdaptiveArithmeticCoding,
            &[0x06],
        )?;
        t(&mut buf, CompressionMethod::Fqzcomp, &[0x07])?;
        t(&mut buf, CompressionMethod::NameTokenizer, &[0x08])?;

        Ok(())
    }

    #[test]
    fn test_write_content_type() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, content_type: ContentType, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_content_type(buf, content_type)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, ContentType::FileHeader, &[0x00])?;
        t(&mut buf, ContentType::CompressionHeader, &[0x01])?;
        t(&mut buf, ContentType::SliceHeader, &[0x02])?;
        t(&mut buf, ContentType::Reserved, &[0x03])?;
        t(&mut buf, ContentType::ExternalData, &[0x04])?;
        t(&mut buf, ContentType::CoreData, &[0x05])?;

        Ok(())
    }
}
