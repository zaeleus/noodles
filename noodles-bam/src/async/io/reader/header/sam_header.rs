//! Async BAM header SAM header reader.

use bstr::ByteSlice;
use pin_project_lite::pin_project;
use std::{
    pin::Pin,
    task::{Context, Poll, ready},
};
use tokio::io::{
    self, AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncReadExt, BufReader, ReadBuf, Take,
};

pin_project! {
    /// An async BAM header SAM header reader.
    pub struct Reader<R> {
        #[pin]
        inner: BufReader<Take<R>>,
        is_eol: bool,
    }
}

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    pub(super) fn new(inner: R, len: u64) -> Self {
        Self {
            inner: BufReader::new(inner.take(len)),
            is_eol: true,
        }
    }

    /// Discards all input until EOF.
    pub async fn discard_to_end(&mut self) -> io::Result<usize> {
        let mut n = 0;

        loop {
            let src = self.inner.fill_buf().await?;

            if src.is_empty() {
                return Ok(n);
            }

            let len = src.len();

            self.inner.consume(len);

            n += len;
        }
    }
}

impl<R> AsyncRead for Reader<R>
where
    R: AsyncRead + Unpin,
{
    fn poll_read(
        mut self: Pin<&mut Self>,
        cx: &mut Context<'_>,
        buf: &mut ReadBuf<'_>,
    ) -> Poll<io::Result<()>> {
        let src = ready!(self.as_mut().poll_fill_buf(cx))?;

        let amt = src.len().min(buf.remaining());
        buf.put_slice(&src[..amt]);

        if amt < src.len() {
            self.is_eol = false;
        }

        self.consume(amt);

        Poll::Ready(Ok(()))
    }
}

impl<R> AsyncBufRead for Reader<R>
where
    R: AsyncRead + Unpin,
{
    fn poll_fill_buf(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<io::Result<&[u8]>> {
        const NUL: u8 = 0x00;
        const LINE_FEED: u8 = b'\n';

        let this = self.project();
        let src = ready!(this.inner.poll_fill_buf(cx))?;

        let buf = if *this.is_eol && src.first().map(|&b| b == NUL).unwrap_or(true) {
            &[]
        } else if let Some(i) = src.as_bstr().find_byte(LINE_FEED) {
            *this.is_eol = true;
            &src[..=i]
        } else {
            *this.is_eol = false;
            src
        };

        Poll::Ready(Ok(buf))
    }

    fn consume(mut self: Pin<&mut Self>, amt: usize) {
        self.as_mut().inner.consume(amt);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_with_trailing_nul_padding() -> io::Result<()> {
        const DATA: &[u8] = b"@HD\tVN:1.6\n";

        let mut buf = DATA.to_vec();
        buf.resize(1 << 10, 0);

        let mut src = &buf[..];
        let mut reader = Reader::new(&mut src, 1 << 10);

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).await?;

        assert_eq!(buf, DATA);

        Ok(())
    }
}
