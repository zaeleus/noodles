use std::io::{self, BufRead};

use crate::Header;

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: BufRead,
{
    read_raw_header(reader).and_then(|s| {
        s.parse()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    })
}

fn read_raw_header<R>(reader: &mut R) -> io::Result<String>
where
    R: BufRead,
{
    use memchr::memchr;

    const HEADER_PREFIX: u8 = b'#';
    const LINE_FEED: u8 = b'\n';

    let mut buf = Vec::new();

    let mut is_first_line = true;
    let mut is_eol = false;

    loop {
        let src = reader.fill_buf()?;

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

    #[test]
    fn test_read_raw_header() -> io::Result<()> {
        static DATA: &[u8] = b"\
##fileformat=VCFv4.3
##fileDate=20200501
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
sq0\t1\t.\tA\t.\t.\tPASS\t.
";

        let mut reader = DATA;

        let actual = read_raw_header(&mut reader)?;
        let expected = "\
##fileformat=VCFv4.3
##fileDate=20200501
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_raw_header_with_no_records() -> io::Result<()> {
        let expected = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

        let mut reader = expected.as_bytes();
        let actual = read_raw_header(&mut reader)?;

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_raw_header_with_multiple_buffer_fills() -> io::Result<()> {
        use std::io::BufReader;

        let expected = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

        let mut reader = BufReader::with_capacity(16, expected.as_bytes());
        let actual = read_raw_header(&mut reader)?;

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_raw_header_with_no_header() -> io::Result<()> {
        let data = [];
        let mut reader = &data[..];
        let actual = read_raw_header(&mut reader)?;
        assert!(actual.is_empty());

        let data = b"sq0\t1\t.\tA\t.\t.\tPASS\t.\n";
        let mut reader = &data[..];
        let actual = read_raw_header(&mut reader)?;
        assert!(actual.is_empty());

        Ok(())
    }

    #[test]
    fn test_read_raw_header_with_missing_end_of_line() -> io::Result<()> {
        let expected = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        let mut reader = expected.as_bytes();
        let actual = read_raw_header(&mut reader)?;
        assert_eq!(actual, expected);
        Ok(())
    }
}
