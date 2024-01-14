mod header;

use std::str;

use bytes::{Buf, Bytes, BytesMut};
use tokio::io::{self, AsyncRead, AsyncReadExt};

use self::header::read_header;
use crate::container::{
    block::{CompressionMethod, ContentType},
    Block,
};

pub async fn read_header_container<R>(reader: &mut R, buf: &mut BytesMut) -> io::Result<String>
where
    R: AsyncRead + Unpin,
{
    let len = read_header(reader).await?;

    buf.resize(len, 0);
    reader.read_exact(buf).await?;
    let mut buf = buf.split().freeze();

    read_raw_sam_header_from_block(&mut buf)
}

fn read_raw_sam_header_from_block(src: &mut Bytes) -> io::Result<String> {
    use crate::io::reader::container::read_block;

    let block = read_block(src)?;
    read_raw_sam_header(&block)
}

fn read_raw_sam_header(block: &Block) -> io::Result<String> {
    const EXPECTED_CONTENT_TYPE: ContentType = ContentType::FileHeader;

    if !matches!(
        block.compression_method(),
        CompressionMethod::None | CompressionMethod::Gzip
    ) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid block compression method: expected {:?} or {:?}, got {:?}",
                CompressionMethod::None,
                CompressionMethod::Gzip,
                block.compression_method()
            ),
        ));
    }

    if block.content_type() != EXPECTED_CONTENT_TYPE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid block content type: expected {:?}, got {:?}",
                EXPECTED_CONTENT_TYPE,
                block.content_type()
            ),
        ));
    }

    let mut data = block.decompressed_data()?;

    let len = usize::try_from(data.get_i32_le())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    data.truncate(len);

    str::from_utf8(&data[..])
        .map(|s| s.into())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use bytes::BufMut;

    use super::*;

    #[test]
    fn test_read_raw_sam_header() -> io::Result<()> {
        let raw_header = "@HD\tVN:1.6\n";

        let header_data = raw_header.to_string().into_bytes();
        let header_data_len = i32::try_from(header_data.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let mut data = Vec::new();
        data.put_i32_le(header_data_len);
        data.extend(&header_data);

        let block = Block::builder()
            .set_content_type(ContentType::FileHeader)
            .set_uncompressed_len(data.len())
            .set_data(data.into())
            .build();

        let actual = read_raw_sam_header(&block)?;

        assert_eq!(actual, raw_header);

        Ok(())
    }

    #[test]
    fn test_read_raw_sam_header_with_invalid_compression_method() {
        let block = Block::builder()
            .set_compression_method(CompressionMethod::Lzma)
            .set_content_type(ContentType::FileHeader)
            .build();

        assert!(matches!(
            read_raw_sam_header(&block),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[test]
    fn test_read_raw_sam_header_with_invalid_content_type() {
        let block = Block::builder()
            .set_content_type(ContentType::ExternalData)
            .build();

        assert!(matches!(
            read_raw_sam_header(&block),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }
}
