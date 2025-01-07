mod block;
mod header;

use bytes::BytesMut;
use tokio::io::{self, AsyncRead, AsyncReadExt};

use self::{block::read_block, header::read_header};

pub async fn read_header_container<R>(reader: &mut R, buf: &mut BytesMut) -> io::Result<String>
where
    R: AsyncRead + Unpin,
{
    let len = read_header(reader).await?;

    buf.resize(len, 0);
    reader.read_exact(buf).await?;

    let buf = buf.split().freeze();
    let mut reader = &buf[..];
    read_raw_sam_header(&mut reader).await
}

async fn read_raw_sam_header<R>(reader: &mut R) -> io::Result<String>
where
    R: AsyncRead + Unpin,
{
    let mut block_reader = read_block(reader).await?;

    let len = block_reader.read_i32_le().await.and_then(|n| {
        u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut header_reader = block_reader.take(len);
    let mut buf = String::new();
    header_reader.read_to_string(&mut buf).await?;

    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_raw_sam_header() -> io::Result<()> {
        const RAW_HEADER: &str = "@HD\tVN:1.6\n";

        let mut src = vec![
            0x00, // compression method = none (0)
            0x00, // content type = file header (0)
            0x00, // block content ID
            0x0f, // compressed size = 15
            0x0f, // uncompressed size = 15
        ];
        src.extend(11i32.to_le_bytes());
        src.extend(RAW_HEADER.as_bytes());
        src.extend([0x00, 0x00, 0x00, 0x00]); // CRC32

        let mut reader = &src[..];
        let actual = read_raw_sam_header(&mut reader).await?;

        assert_eq!(actual, RAW_HEADER);

        Ok(())
    }

    #[tokio::test]
    async fn test_read_raw_sam_header_with_invalid_compression_method() {
        let src = [
            0x03, // compression method = LZMA (3)
            0x00, // content type = file header (0)
            0x00, // block content ID
            0x0f, // compressed size = 15
            0x0f, // uncompressed size = 15
            // ...
            0x00, 0x00, 0x00, 0x00, // CRC32
        ];

        let mut reader = &src[..];
        assert!(matches!(
            read_raw_sam_header(&mut reader).await,
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[tokio::test]
    async fn test_read_raw_sam_header_with_invalid_content_type() {
        let src = [
            0x00, // compression method = none (0)
            0x04, // content type = external data (4)
        ];

        let mut reader = &src[..];

        assert!(matches!(
            read_raw_sam_header(&mut reader).await,
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }
}
