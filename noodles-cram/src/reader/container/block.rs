use std::{io, mem};

use bytes::{Buf, Bytes};

use crate::{
    container::{
        block::{CompressionMethod, ContentId, ContentType},
        Block,
    },
    reader::num::get_itf8,
};

pub fn read_block(src: &mut Bytes) -> io::Result<Block> {
    let original_src = src.clone();

    let method = get_compression_method(src)?;

    if !src.has_remaining() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let block_content_type_id = ContentType::try_from(src.get_u8())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let block_content_id = get_itf8(src).map(ContentId::from)?;

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
        .build())
}

fn get_compression_method<B>(src: &mut B) -> io::Result<CompressionMethod>
where
    B: Buf,
{
    if !src.has_remaining() {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    match src.get_u8() {
        0 => Ok(CompressionMethod::None),
        1 => Ok(CompressionMethod::Gzip),
        2 => Ok(CompressionMethod::Bzip2),
        3 => Ok(CompressionMethod::Lzma),
        4 => Ok(CompressionMethod::Rans4x8),
        5 => Ok(CompressionMethod::RansNx16),
        6 => Ok(CompressionMethod::AdaptiveArithmeticCoding),
        7 => Ok(CompressionMethod::Fqzcomp),
        8 => Ok(CompressionMethod::NameTokenizer),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid compression method",
        )),
    }
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
            .set_content_id(ContentId::from(1))
            .set_uncompressed_len(4)
            .set_data(Bytes::from_static(b"ndls"))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_get_compression_method() -> io::Result<()> {
        fn t(mut src: &[u8], expected: CompressionMethod) -> io::Result<()> {
            let actual = get_compression_method(&mut src)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(&[0x00], CompressionMethod::None)?;
        t(&[0x01], CompressionMethod::Gzip)?;
        t(&[0x02], CompressionMethod::Bzip2)?;
        t(&[0x03], CompressionMethod::Lzma)?;
        t(&[0x04], CompressionMethod::Rans4x8)?;
        t(&[0x05], CompressionMethod::RansNx16)?;
        t(&[0x06], CompressionMethod::AdaptiveArithmeticCoding)?;
        t(&[0x07], CompressionMethod::Fqzcomp)?;
        t(&[0x08], CompressionMethod::NameTokenizer)?;

        let mut src = &[][..];
        assert!(matches!(
            get_compression_method(&mut src),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let mut src = &[0x09][..];
        assert!(matches!(
            get_compression_method(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
