//! BAM header SAM header reader.

use std::io::{self, BufRead, BufReader, Read, Take};

use bstr::ByteSlice;

/// A BAM header SAM header reader.
pub struct Reader<R> {
    inner: BufReader<Take<R>>,
    is_eol: bool,
}

impl<R> Reader<R>
where
    R: Read,
{
    pub(super) fn new(inner: R, len: u64) -> Self {
        Self {
            inner: BufReader::new(inner.take(len)),
            is_eol: true,
        }
    }

    /// Discards all input until EOF.
    pub fn discard_to_end(&mut self) -> io::Result<usize> {
        let mut n = 0;

        loop {
            let src = self.inner.fill_buf()?;

            if src.is_empty() {
                return Ok(n);
            }

            let len = src.len();

            self.inner.consume(len);

            n += len;
        }
    }
}

impl<R> Read for Reader<R>
where
    R: Read,
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
    R: Read,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        const NUL: u8 = 0x00;
        const LINE_FEED: u8 = b'\n';

        let src = self.inner.fill_buf()?;

        let buf = if self.is_eol && src.first().map(|&b| b == NUL).unwrap_or(true) {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_with_trailing_nul_padding() -> io::Result<()> {
        const DATA: &[u8] = b"@HD\tVN:1.6\n";

        let mut buf = DATA.to_vec();
        buf.resize(1 << 10, 0);

        let mut src = &buf[..];
        let mut reader = Reader::new(&mut src, 1 << 10);

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf)?;

        assert_eq!(buf, DATA);

        Ok(())
    }
}
