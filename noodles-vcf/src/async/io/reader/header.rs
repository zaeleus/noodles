use std::{
    pin::Pin,
    task::{ready, Context, Poll},
};

use pin_project_lite::pin_project;
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncReadExt, ReadBuf};

use crate::Header;

pin_project! {
    struct Reader<R> {
        inner: R,
        is_eol: bool,
    }
}

impl<R> Reader<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            is_eol: true,
        }
    }
}

impl<R> AsyncRead for Reader<R>
where
    R: AsyncBufRead + Unpin,
{
    fn poll_read(
        mut self: Pin<&mut Self>,
        cx: &mut Context<'_>,
        buf: &mut ReadBuf<'_>,
    ) -> Poll<io::Result<()>> {
        let src = ready!(self.as_mut().poll_fill_buf(cx))?;

        let amt = src.len().min(buf.remaining());
        buf.put_slice(&src[..amt]);

        self.consume(amt);

        Poll::Ready(Ok(()))
    }
}

impl<R> AsyncBufRead for Reader<R>
where
    R: AsyncBufRead + Unpin,
{
    fn poll_fill_buf(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<io::Result<&[u8]>> {
        use memchr::memchr;

        const PREFIX: u8 = b'#';
        const LINE_FEED: u8 = b'\n';

        let this = self.project();
        let buf = ready!(Pin::new(this.inner).poll_fill_buf(cx))?;

        let buf = if *this.is_eol && buf.first().map(|&b| b != PREFIX).unwrap_or(true) {
            &[]
        } else if let Some(i) = memchr(LINE_FEED, buf) {
            *this.is_eol = true;
            &buf[..=i]
        } else {
            *this.is_eol = false;
            buf
        };

        Poll::Ready(Ok(buf))
    }

    fn consume(self: Pin<&mut Self>, amt: usize) {
        let this = self.project();
        this.inner.consume(amt);
    }
}

pub(super) async fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: AsyncBufRead + Unpin,
{
    let mut s = String::new();
    Reader::new(reader).read_to_string(&mut s).await?;
    s.parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
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

        let mut src = DATA;
        let mut reader = Reader::new(&mut src);

        let mut actual = String::new();
        reader.read_to_string(&mut actual).await?;

        let expected = "\
##fileformat=VCFv4.3
##fileDate=20200501
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        assert_eq!(actual, expected);

        Ok(())
    }
}
