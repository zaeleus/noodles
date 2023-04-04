use tokio::io::{self, AsyncBufRead, AsyncBufReadExt};

use crate::Header;

pub(super) async fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: AsyncBufRead + Unpin,
{
    read_raw_header(reader).await.and_then(|s| {
        s.parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    })
}

async fn read_raw_header<R>(reader: &mut R) -> io::Result<String>
where
    R: AsyncBufRead + Unpin,
{
    use memchr::memchr;

    const HEADER_PREFIX: u8 = b'#';
    const LINE_FEED: u8 = b'\n';

    let mut buf = Vec::new();

    let mut is_first_line = true;
    let mut is_eol = false;

    loop {
        let src = reader.fill_buf().await?;

        let is_eof = src.is_empty();
        let is_end_of_header = || (is_first_line || is_eol) && src[0] != HEADER_PREFIX;

        if is_eof || is_end_of_header() {
            break;
        }

        let (read_eol, len) = if let Some(i) = memchr(LINE_FEED, src) {
            buf.extend(&src[..=i]);
            (true, i + 1)
        } else {
            buf.extend(src);
            (false, src.len())
        };

        is_first_line = false;
        is_eol = read_eol;

        reader.consume(len);
    }

    String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_raw_header() -> io::Result<()> {
        static DATA: &[u8] = b"\
##fileformat=VCFv4.3
##fileDate=20200501
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
sq0\t1\t.\tA\t.\t.\tPASS\t.
";

        let mut reader = DATA;

        let actual = read_raw_header(&mut reader).await?;
        let expected = "\
##fileformat=VCFv4.3
##fileDate=20200501
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        assert_eq!(actual, expected);

        Ok(())
    }
}
