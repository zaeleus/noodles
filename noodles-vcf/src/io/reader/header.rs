//! VCF header reader.

use std::io::{self, BufRead, Read};

use crate::{header, Header};

/// A VCF header reader.
///
/// This is created by calling [`super::Reader::header_reader`].
pub struct Reader<R> {
    inner: R,
    is_eol: bool,
}

impl<R> Reader<R> {
    pub(super) fn new(inner: R) -> Self {
        Self {
            inner,
            is_eol: true,
        }
    }
}

impl<R> Read for Reader<R>
where
    R: BufRead,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut src = self.fill_buf()?;
        let amt = src.read(buf)?;

        if !src.is_empty() {
            self.is_eol = false;
        }

        self.consume(amt);

        Ok(amt)
    }
}

impl<R> BufRead for Reader<R>
where
    R: BufRead,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        use memchr::memchr;

        const PREFIX: u8 = b'#';
        const LINE_FEED: u8 = b'\n';

        let buf = self.inner.fill_buf()?;

        if self.is_eol && buf.first().map(|&b| b != PREFIX).unwrap_or(true) {
            Ok(&[])
        } else if let Some(i) = memchr(LINE_FEED, buf) {
            self.is_eol = true;
            Ok(&buf[..=i])
        } else {
            self.is_eol = false;
            Ok(buf)
        }
    }

    fn consume(&mut self, amt: usize) {
        self.inner.consume(amt);
    }
}

pub(super) fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: BufRead,
{
    let mut reader = Reader::new(reader);

    let mut parser = header::Parser::default();
    let mut buf = Vec::new();

    while read_line(&mut reader, &mut buf)? != 0 {
        parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    parser
        .finish()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_line<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    dst.clear();

    match reader.read_until(LINE_FEED, dst)? {
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

    #[test]
    fn test_read_raw_header() -> io::Result<()> {
        let data = b"##fileformat=VCFv4.3
##fileDate=20200501
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
sq0\t1\t.\tA\t.\t.\tPASS\t.
";

        let mut src = &data[..];
        let mut reader = Reader::new(&mut src);

        let mut actual = Vec::new();
        reader.read_to_end(&mut actual)?;

        let expected = b"##fileformat=VCFv4.3
##fileDate=20200501
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_raw_header_with_no_records() -> io::Result<()> {
        let data = b"##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        let mut src = &data[..];
        let mut reader = Reader::new(&mut src);

        let mut actual = Vec::new();
        reader.read_to_end(&mut actual)?;

        let expected = b"##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_raw_header_with_multiple_buffer_fills() -> io::Result<()> {
        use std::io::BufReader;

        let data = b"##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        let mut inner = BufReader::with_capacity(16, &data[..]);
        let mut reader = Reader::new(&mut inner);

        let mut actual = Vec::new();
        reader.read_to_end(&mut actual)?;

        let expected = b"##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_raw_header_with_no_header() -> io::Result<()> {
        let mut src = &[][..];
        let mut reader = Reader::new(&mut src);
        let mut actual = Vec::new();
        reader.read_to_end(&mut actual)?;
        assert!(actual.is_empty());

        let mut src = &b"sq0\t1\t.\tA\t.\t.\tPASS\t.\n"[..];
        let mut reader = Reader::new(&mut src);
        let mut actual = Vec::new();
        reader.read_to_end(&mut actual)?;
        assert!(actual.is_empty());

        Ok(())
    }

    #[test]
    fn test_read_raw_header_with_missing_end_of_line() -> io::Result<()> {
        let data = b"##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

        let mut src = &data[..];
        let mut reader = Reader::new(&mut src);

        let mut actual = Vec::new();
        reader.read_to_end(&mut actual)?;

        let expected = b"##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

        assert_eq!(actual, expected);

        Ok(())
    }
}
