use std::convert::TryFrom;

use tokio::io::{self, AsyncRead, AsyncReadExt};

use super::num::read_itf8;
use crate::container::{
    block::{CompressionMethod, ContentType},
    Block,
};

pub async fn read_block<R>(reader: &mut R) -> io::Result<Block>
where
    R: AsyncRead + Unpin,
{
    let method = reader.read_u8().await.and_then(|n| {
        CompressionMethod::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let block_content_type_id = reader.read_u8().await.and_then(|n| {
        ContentType::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let block_content_id = read_itf8(reader).await?;
    let size_in_bytes = read_itf8(reader).await?;
    let raw_size_in_bytes = read_itf8(reader).await?;

    let len = usize::try_from(size_in_bytes)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut block_data = vec![0; len];
    reader.read_exact(&mut block_data).await?;

    let crc32 = reader.read_u32_le().await?;

    Ok(Block::builder()
        .set_compression_method(method)
        .set_content_type(block_content_type_id)
        .set_content_id(block_content_id)
        .set_uncompressed_len(raw_size_in_bytes)
        .set_data(block_data)
        .set_crc32(crc32)
        .build())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_block() -> io::Result<()> {
        let data = [
            0x00, // compression method = none (0)
            0x04, // content type = EXTERNAL_DATA (4)
            0x01, // block content ID = 1
            0x04, // size in bytes = 4 bytes
            0x04, // raw size in bytes = 4 bytes
            0x6e, 0x64, 0x6c, 0x73, // block data = b"ndls",
            0xfd, 0x38, 0x27, 0xb5, // CRC32
        ];

        let mut reader = &data[..];
        let actual = read_block(&mut reader).await?;

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
