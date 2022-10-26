mod header;

use tokio::io::{self, AsyncWrite};

use bytes::BufMut;
use noodles_sam as sam;

use self::header::write_header;
use super::container::write_block;
use crate::container::{block::ContentType, Block};

pub async fn write_header_container<W>(writer: &mut W, header: &sam::Header) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let header_data = header.to_string().into_bytes();
    let header_data_len = i32::try_from(header_data.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let mut data = Vec::new();
    data.put_i32_le(header_data_len);
    data.extend_from_slice(&header_data);

    let block = Block::builder()
        .set_content_type(ContentType::FileHeader)
        .set_uncompressed_len(data.len())
        .set_data(data.into())
        .build();

    write_header(writer, block.len()).await?;
    write_block(writer, &block).await?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_header_container() -> io::Result<()> {
        let header = sam::Header::builder()
            .set_header(Default::default())
            .build();

        let mut buf = Vec::new();
        write_header_container(&mut buf, &header).await?;

        let mut expected = Vec::new();

        // header
        expected.extend_from_slice(&[
            0x18, 0x00, 0x00, 0x00, // length = 24
            0xff, 0xff, 0xff, 0xff, 0x0f, // reference sequence ID = -1 (None)
            0x00, // alignment start = 0
            0x00, // alignment span = 0
            0x00, // record count = 0
            0x00, // record counter = 0
            0x00, // base count = 0
            0x01, // block count = 1
            0x00, // landmarks.len = 0
            0xf9, 0xf0, 0x9e, 0x3d, // CRC32 = 3d9ef0f9
        ]);

        // block
        expected.extend_from_slice(&[
            0x00, // compression method = 0 (None)
            0x00, // content type = 0 (FileHeader)
            0x00, // block content ID = 0
            0x0f, // compressed length = 15
            0x0f, // uncompressed length = 15
        ]);
        expected.put_i32_le(11);
        expected.extend_from_slice(b"@HD\tVN:1.6\n");
        expected.extend_from_slice(&[
            0xf4, 0xe5, 0x16, 0xd3, // CRC32 = d316e5f4
        ]);

        assert_eq!(buf, expected);

        Ok(())
    }
}
