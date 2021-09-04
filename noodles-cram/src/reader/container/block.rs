use std::{
    convert::TryFrom,
    io::{self, Read},
};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::{
    container::{
        block::{CompressionMethod, ContentType},
        Block,
    },
    reader::num::read_itf8,
};

pub fn read_block<R>(reader: &mut R) -> io::Result<Block>
where
    R: Read,
{
    let method = reader.read_u8().and_then(|b| {
        CompressionMethod::try_from(b).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let block_content_type_id = reader.read_u8().and_then(|b| {
        ContentType::try_from(b).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let block_content_id = read_itf8(reader)?;

    let size_in_bytes = read_itf8(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let raw_size_in_bytes = read_itf8(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut data = vec![0; size_in_bytes];
    reader.read_exact(&mut data)?;

    let crc32 = reader.read_u32::<LittleEndian>()?;

    Ok(Block::builder()
        .set_compression_method(method)
        .set_content_type(block_content_type_id)
        .set_content_id(block_content_id)
        .set_uncompressed_len(raw_size_in_bytes)
        .set_data(data)
        .set_crc32(crc32)
        .build())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_block() -> io::Result<()> {
        let data = [
            0x00, // compression method = none (0)
            0x04, // content type = external data (4)
            0x01, // block content ID = 1
            0x04, // size in bytes = 4 bytes
            0x04, // raw size in bytes = 4 bytes
            0x6e, 0x64, 0x6c, 0x73, // data = b"ndls",
            0xfd, 0x38, 0x27, 0xb5, // CRC32
        ];
        let mut reader = &data[..];
        let actual = read_block(&mut reader)?;

        let expected = Block::builder()
            .set_compression_method(CompressionMethod::None)
            .set_content_type(ContentType::ExternalData)
            .set_content_id(1)
            .set_uncompressed_len(4)
            .set_data(b"ndls".to_vec())
            .set_crc32(0xb52738fd)
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
