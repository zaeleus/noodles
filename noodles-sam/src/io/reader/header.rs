//! SAM header reader.

use std::io::{self, BufRead, Read};

use bstr::ByteSlice;

use super::read_line;
use crate::{Header, header};

/// A SAM header reader.
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
        const PREFIX: u8 = b'@';
        const LINE_FEED: u8 = b'\n';

        let src = self.inner.fill_buf()?;

        let buf = if self.is_eol && src.first().map(|&b| b != PREFIX).unwrap_or(true) {
            &[]
        } else if let Some(i) = src.as_bstr().find_byte(LINE_FEED) {
            self.is_eol = true;
            &src[..=i]
        } else {
            self.is_eol = false;
            src
        };

        Ok(buf)
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

        buf.clear();
    }

    Ok(parser.finish())
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use super::*;
    use crate::header::record::value::{
        Map,
        map::{self, ReferenceSequence, header::Version},
    };

    #[test]
    fn test_read_header_with_no_header() -> io::Result<()> {
        let data = b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n";
        let mut reader = &data[..];
        assert!(read_header(&mut reader)?.is_empty());
        Ok(())
    }

    #[test]
    fn test_read_header_with_no_records() -> io::Result<()> {
        let data = "@HD\tVN:1.6\n";
        let mut reader = data.as_bytes();

        let actual = read_header(&mut reader)?;

        let expected = crate::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_header_with_multiple_buffer_fills() -> io::Result<()> {
        use std::io::BufReader;

        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        let data = "@HD\tVN:1.6\n@SQ\tSN:sq0\tLN:8\n";
        let mut reader = BufReader::with_capacity(16, data.as_bytes());

        let actual = read_header(&mut reader)?;

        let expected = crate::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(SQ0_LN))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
