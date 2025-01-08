mod block;
mod header;

use noodles_sam as sam;
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncReadExt, BufReader};

use self::{block::read_block, header::read_header};

pub async fn read_header_container<R>(reader: &mut R) -> io::Result<sam::Header>
where
    R: AsyncRead + Unpin,
{
    let len = read_header(reader).await?;

    let mut reader = reader.take(len);
    let header = read_sam_header(&mut reader).await?;
    io::copy(&mut reader, &mut io::sink()).await?;

    Ok(header)
}

async fn read_sam_header<R>(reader: &mut R) -> io::Result<sam::Header>
where
    R: AsyncRead + Unpin,
{
    let mut block_reader = read_block(reader).await?;

    let len = block_reader.read_i32_le().await.and_then(|n| {
        u64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut parser = sam::header::Parser::default();

    let mut header_reader = BufReader::new(block_reader.take(len));
    let mut buf = Vec::new();

    while read_line(&mut header_reader, &mut buf).await? != 0 {
        parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    Ok(parser.finish())
}

async fn read_line<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    dst.clear();

    match reader.read_until(LINE_FEED, dst).await? {
        0 => Ok(0),
        n => {
            if dst.ends_with(&[LINE_FEED]) {
                dst.pop();

                if dst.ends_with(&[CARRIAGE_RETURN]) {
                    dst.pop();
                }
            }

            Ok(n)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_raw_sam_header() -> io::Result<()> {
        use sam::header::record::value::{
            map::{self, header::Version},
            Map,
        };

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
        let actual = read_sam_header(&mut reader).await?;

        let expected = sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .build();

        assert_eq!(actual, expected);

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
            read_sam_header(&mut reader).await,
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
            read_sam_header(&mut reader).await,
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }
}
