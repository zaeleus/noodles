//! Async SAM header reader.

use bstr::ByteSlice;
use pin_project_lite::pin_project;
use std::{
    pin::Pin,
    task::{ready, Context, Poll},
};
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, ReadBuf};

use super::read_line;
use crate::{header, Header};

pin_project! {
    /// An async SAM header reader.
    ///
    /// This is created by calling [`super::Reader::header_reader`].
    pub struct Reader<R> {
        #[pin]
        inner: R,
        is_eol: bool,
    }
}

impl<R> Reader<R> {
    pub(super) fn new(inner: R) -> Self {
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
        const PREFIX: u8 = b'@';
        const LINE_FEED: u8 = b'\n';

        let this = self.project();
        let src = ready!(this.inner.poll_fill_buf(cx))?;

        let buf = if *this.is_eol && src.first().map(|&b| b != PREFIX).unwrap_or(true) {
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

        buf.clear();
    }

    Ok(parser.finish())
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use super::*;
    use crate::header::record::value::{
        map::{self, header::Version, ReferenceSequence},
        Map,
    };

    #[tokio::test]
    async fn test_read_header_with_no_header() -> io::Result<()> {
        let data = b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n";
        let mut reader = &data[..];
        assert!(read_header(&mut reader).await?.is_empty());
        Ok(())
    }

    #[tokio::test]
    async fn test_read_header_with_no_records() -> io::Result<()> {
        let data = "@HD\tVN:1.6\n";
        let mut reader = data.as_bytes();

        let actual = read_header(&mut reader).await?;

        let expected = crate::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[tokio::test]
    async fn test_read_header_with_multiple_buffer_fills() -> io::Result<()> {
        use tokio::io::BufReader;

        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        let data = "@HD\tVN:1.6\n@SQ\tSN:sq0\tLN:8\n";
        let mut reader = BufReader::with_capacity(16, data.as_bytes());

        let actual = read_header(&mut reader).await?;

        let expected = crate::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(SQ0_LN))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
