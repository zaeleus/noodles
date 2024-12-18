use std::io::{self, BufRead, BufReader, Read, Take};

use bstr::ByteSlice;

pub(super) struct Reader<'r, R> {
    inner: BufReader<Take<&'r mut R>>,
    is_eol: bool,
}

impl<'r, R> Reader<'r, R>
where
    R: Read,
{
    pub(super) fn new(inner: &'r mut R, len: u64) -> Self {
        Self {
            inner: BufReader::new(inner.take(len)),
            is_eol: true,
        }
    }

    pub(super) fn discard_to_end(&mut self) -> io::Result<usize> {
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

impl<R> Read for Reader<'_, R>
where
    R: Read,
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
    R: Read,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        const NUL: u8 = 0x00;
        const LINE_FEED: u8 = b'\n';

        let src = self.inner.fill_buf()?;

        if self.is_eol && src.first().map(|&b| b == NUL).unwrap_or(true) {
            Ok(&[])
        } else if let Some(i) = src.as_bstr().find_byte(LINE_FEED) {
            self.is_eol = true;
            Ok(&src[..=i])
        } else {
            self.is_eol = false;
            Ok(src)
        }
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
