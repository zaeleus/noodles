//! Sequence reader.

use std::io::{self, BufRead, Read};

use super::DEFINITION_PREFIX;

const LINE_FEED: u8 = b'\n';
const CARRIAGE_RETURN: u8 = b'\r';

/// A sequence reader.
///
/// This is used for lower-level reading of the sequence. It implements [`Read`] and [`BufRead`] to
/// return raw bases sans newlines. It reads up to the next record definition or EOF.
///
/// This is created by calling [`super::Reader::sequence_reader`].
pub struct Reader<'r, R> {
    inner: &'r mut R,
}

impl<'r, R> Reader<'r, R>
where
    R: BufRead,
{
    pub(super) fn new(inner: &'r mut R) -> Self {
        Self { inner }
    }
}

impl<R> Read for Reader<'_, R>
where
    R: BufRead,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut src = self.fill_buf()?;
        let amt = src.read(buf)?;
        self.consume(amt);
        Ok(amt)
    }
}

impl<R> BufRead for Reader<'_, R>
where
    R: BufRead,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        use memchr::memchr;

        consume_empty_lines(&mut self.inner)?;

        let src = self.inner.fill_buf()?;

        let is_eof = src.is_empty();
        let is_end_of_sequence = || src[0] == DEFINITION_PREFIX;

        if is_eof || is_end_of_sequence() {
            return Ok(&[]);
        }

        let line = match memchr(LINE_FEED, src) {
            Some(i) => &src[..i],
            None => src,
        };

        if line.ends_with(&[CARRIAGE_RETURN]) {
            let end = line.len() - 1;
            Ok(&line[..end])
        } else {
            Ok(line)
        }
    }

    fn consume(&mut self, amt: usize) {
        self.inner.consume(amt);
    }
}

fn consume_empty_lines<R>(reader: &mut R) -> io::Result<()>
where
    R: BufRead,
{
    loop {
        let mut is_newline = false;

        if reader.fill_buf()?.starts_with(&[CARRIAGE_RETURN]) {
            is_newline = true;
            reader.consume(1);
        }

        if reader.fill_buf()?.starts_with(&[LINE_FEED]) {
            is_newline = true;
            reader.consume(1);
        }

        if !is_newline {
            break;
        }
    }

    Ok(())
}

pub(super) fn read_sequence<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: BufRead,
{
    let mut reader = Reader::new(reader);
    reader.read_to_end(buf)
}

pub(super) fn read_sequence_limit<R>(
    reader: &mut R,
    max_bases: usize,
    buf: &mut Vec<u8>,
) -> io::Result<usize>
where
    R: BufRead,
{
    let mut reader = Reader::new(reader);
    let mut len = 0;

    while buf.len() < max_bases {
        let src = reader.fill_buf()?;

        if src.is_empty() {
            break;
        }

        let remaining_bases = max_bases - buf.len();
        let i = remaining_bases.min(src.len());

        let bases = &src[..i];
        buf.extend(bases);

        reader.consume(i);

        len += i;
    }

    Ok(len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_sequence() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, mut reader: &[u8], expected: &[u8]) -> io::Result<()> {
            buf.clear();
            read_sequence(&mut reader, buf)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, b"ACGT\n", b"ACGT")?;
        t(&mut buf, b"ACGT\n>sq1\n", b"ACGT")?;
        t(&mut buf, b"ACGT\n\nACGT\nAC\n\n", b"ACGTACGTAC")?;

        t(&mut buf, b"ACGT\r\n", b"ACGT")?;
        t(&mut buf, b"ACGT\r\n>sq1\r\n", b"ACGT")?;
        t(&mut buf, b"ACGT\r\n\r\nACGT\r\nAC\r\n\r\n", b"ACGTACGTAC")?;

        Ok(())
    }

    #[test]
    fn test_read_sequence_limit() -> io::Result<()> {
        fn t(
            buf: &mut Vec<u8>,
            mut reader: &[u8],
            max_bases: usize,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            read_sequence_limit(&mut reader, max_bases, buf)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, b"ACGT\n", 4, b"ACGT")?;
        t(&mut buf, b"ACGT\n>sq0\n", 4, b"ACGT")?;
        t(&mut buf, b"ACGT\nACGT\nAC\n", 10, b"ACGTACGTAC")?;

        t(&mut buf, b"ACGT\n", 2, b"AC")?;
        t(&mut buf, b"ACGT\n>sq0\n", 2, b"AC")?;
        t(&mut buf, b"ACGT\nACGT\nAC", 2, b"AC")?;

        t(&mut buf, b"ACGT\n", 5, b"ACGT")?;
        t(&mut buf, b"ACGT\n>sq0\n", 5, b"ACGT")?;
        t(&mut buf, b"ACGT\nACGT\nAC", 5, b"ACGTA")?;

        t(&mut buf, b"ACGT\n", 5, b"ACGT")?;
        t(&mut buf, b"ACGT\r\n>sq0\r\n", 5, b"ACGT")?;
        t(&mut buf, b"ACGT\r\nACGT\r\nAC\r\n", 5, b"ACGTA")?;

        Ok(())
    }
}
