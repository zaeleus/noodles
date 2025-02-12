mod compression_method;

use std::{io, mem};

use bytes::{Buf, Bytes};

use self::compression_method::get_compression_method;
use crate::{
    container::{block::ContentType, Block},
    io::reader::num::get_itf8,
};

pub fn read_block(src: &mut Bytes) -> io::Result<Block> {
    let original_src = src.clone();

    let method = get_compression_method(src)?;
    let block_content_type = get_content_type(src)?;
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
                "container block checksum mismatch: expected {expected_crc32:08x}, got {actual_crc32:08x}"
            ),
        ));
    }

    let mut builder = Block::builder()
        .set_content_type(block_content_type)
        .set_content_id(block_content_id);

    if raw_size_in_bytes > 0 {
        builder = builder
            .set_compression_method(method)
            .set_uncompressed_len(raw_size_in_bytes)
            .set_data(data);
    }

    Ok(builder.build())
}

fn get_content_type<B>(src: &mut B) -> io::Result<ContentType>
where
    B: Buf,
{
    let n = src
        .try_get_u8()
        .map_err(|e| io::Error::new(io::ErrorKind::UnexpectedEof, e))?;

    match n {
        0 => Ok(ContentType::FileHeader),
        1 => Ok(ContentType::CompressionHeader),
        2 => Ok(ContentType::SliceHeader),
        3 => Ok(ContentType::Reserved),
        4 => Ok(ContentType::ExternalData),
        5 => Ok(ContentType::CoreData),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid content type",
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
    use crate::container::block::CompressionMethod;

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
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_block_with_empty_block() -> io::Result<()> {
        let mut data = Bytes::from_static(&[
            0x04, // compression method = rANS 4x8 (4)
            0x04, // content type = external data (4)
            0x01, // block content ID = 1
            0x00, // size in bytes = 0 bytes
            0x00, // raw size in bytes = 0 bytes
            // data = b"",
            0xbd, 0xac, 0x02, 0xbd, // CRC32 = bd02acbd
        ]);
        let actual = read_block(&mut data)?;

        let expected = Block::builder()
            .set_content_type(ContentType::ExternalData)
            .set_content_id(1)
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_get_content_type() -> io::Result<()> {
        fn t(mut src: &[u8], expected: ContentType) -> io::Result<()> {
            let actual = get_content_type(&mut src)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(&[0x00], ContentType::FileHeader)?;
        t(&[0x01], ContentType::CompressionHeader)?;
        t(&[0x02], ContentType::SliceHeader)?;
        t(&[0x03], ContentType::Reserved)?;
        t(&[0x04], ContentType::ExternalData)?;
        t(&[0x05], ContentType::CoreData)?;

        let mut src = &[][..];
        assert!(matches!(
            get_content_type(&mut src),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let mut src = &[0x06][..];
        assert!(matches!(
            get_content_type(&mut src),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
