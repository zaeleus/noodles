use tokio::io::{self, AsyncBufRead, AsyncBufReadExt};

use super::read_line;
use crate::{header, Header};

pub(super) async fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: AsyncBufRead + Unpin,
{
    let mut parser = header::Parser::default();
    let mut buf = Vec::new();

    while read_header_line(reader, &mut buf).await? != 0 {
        parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    Ok(parser.finish())
}

async fn read_header_line<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    const PREFIX: u8 = b'@';

    let src = reader.fill_buf().await?;

    if src.is_empty() || src[0] != PREFIX {
        return Ok(0);
    }

    dst.clear();
    read_line(reader, dst).await
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use super::*;
    use crate::header::record::value::{
        map::{self, header::Version, ReferenceSequence},
        Map,
    };

    #[tokio::test]
    async fn test_read_header_with_no_header() -> io::Result<()> {
        let data = b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n";
        let mut reader = &data[..];
        assert!(read_header(&mut reader).await?.is_empty());
        Ok(())
    }

    #[tokio::test]
    async fn test_read_header_with_no_records() -> io::Result<()> {
        let data = "@HD\tVN:1.6\n";
        let mut reader = data.as_bytes();

        let actual = read_header(&mut reader).await?;

        let expected = crate::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[tokio::test]
    async fn test_read_header_with_multiple_buffer_fills() -> io::Result<()> {
        use tokio::io::BufReader;

        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        let data = "@HD\tVN:1.6\n@SQ\tSN:sq0\tLN:8\n";
        let mut reader = BufReader::with_capacity(16, data.as_bytes());

        let actual = read_header(&mut reader).await?;

        let expected = crate::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(SQ0_LN))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
