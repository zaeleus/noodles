use std::{io, mem};

use bytes::{Buf, Bytes};

use crate::{
    container::{
        block::{CompressionMethod, ContentType},
        Block,
    },
    reader::num::get_itf8,
};

pub fn read_block(src: &mut Bytes) -> io::Result<Block> {
    let original_src = src.clone();

    if !src.has_remaining() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let method = CompressionMethod::try_from(src.get_u8())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    if !src.has_remaining() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let block_content_type_id = ContentType::try_from(src.get_u8())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let block_content_id = get_itf8(src)?;

    let size_in_bytes = get_itf8(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let raw_size_in_bytes = get_itf8(src).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if src.remaining() < size_in_bytes {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let data = src.split_to(size_in_bytes);

    if src.remaining() < mem::size_of::<u32>() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let end = original_src.len() - src.len();
    let actual_crc32 = crc32(&original_src[..end]);

    let expected_crc32 = src.get_u32_le();

    if actual_crc32 != expected_crc32 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "container block checksum mismatch: expected {:08x}, got {:08x}",
                expected_crc32, actual_crc32
            ),
        ));
    }

    Ok(Block::builder()
        .set_compression_method(method)
        .set_content_type(block_content_type_id)
        .set_content_id(block_content_id)
        .set_uncompressed_len(raw_size_in_bytes)
        .set_data(data)
        .set_crc32(expected_crc32)
        .build())
}

fn crc32(buf: &[u8]) -> u32 {
    use flate2::Crc;

    let mut crc = Crc::new();
    crc.update(buf);
    crc.sum()
}

#[cfg(test)]
mod tests {
    use bytes::Bytes;

    use super::*;

    #[test]
    fn test_read_block() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            0x00, // compression method = none (0)
            0x04, // content type = external data (4)
            0x01, // block content ID = 1
            0x04, // size in bytes = 4 bytes
            0x04, // raw size in bytes = 4 bytes
            0x6e, 0x64, 0x6c, 0x73, // data = b"ndls",
            0xd7, 0x12, 0x46, 0x3e, // CRC32 = 3e4612d7
        ]);
        let actual = read_block(&mut data)?;

        let expected = Block::builder()
            .set_compression_method(CompressionMethod::None)
            .set_content_type(ContentType::ExternalData)
            .set_content_id(1)
            .set_uncompressed_len(4)
            .set_data(Bytes::from_static(b"ndls"))
            .set_crc32(0x3e4612d7)
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
