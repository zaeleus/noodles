use std::{
    pin::Pin,
    task::{ready, Context, Poll},
};

use pin_project_lite::pin_project;
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, ReadBuf};

use crate::{header, Header};

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

        if amt < src.len() {
            self.is_eol = false;
        }

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
    let mut reader = Reader::new(reader);

    let mut parser = header::Parser::default();
    let mut buf = Vec::new();

    while read_line(&mut reader, &mut buf).await? != 0 {
        parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    parser
        .finish()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

async fn read_line<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    dst.clear();

    match reader.read_until(LINE_FEED, dst).await? {
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
    use tokio::io::AsyncReadExt;

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
