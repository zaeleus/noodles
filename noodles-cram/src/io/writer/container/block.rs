mod compression_method;
mod content_type;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use flate2::CrcWriter;

use self::{compression_method::write_compression_method, content_type::write_content_type};
use crate::{container::Block, io::writer::num::write_itf8};

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
    write_compression_method(writer, block.compression_method())?;

    write_content_type(writer, block.content_type())?;
    write_itf8(writer, block.content_id())?;

    write_size(writer, block.data().len())?; // compressed size
    write_size(writer, block.uncompressed_len())?;

    writer.write_all(block.data())?;

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

#[cfg(test)]
mod tests {
    use bytes::Bytes;

    use super::*;
    use crate::container::block::{CompressionMethod, ContentType};

    #[test]
    fn test_write_block() -> io::Result<()> {
        let data = Bytes::from_static(b"ndls");

        let block = Block::builder()
            .set_compression_method(CompressionMethod::None)
            .set_content_type(ContentType::ExternalData)
            .set_content_id(1)
            .set_uncompressed_len(data.len())
            .set_data(data)
            .build();

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
}
