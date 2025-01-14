//! Async BGZF writer.

mod builder;
pub(crate) mod deflate;
mod deflater;

use std::{
    pin::Pin,
    task::{ready, Context, Poll},
};

use bytes::{Buf, Bytes, BytesMut};
use futures::{sink::Buffer, Sink};
use pin_project_lite::pin_project;
use tokio::io::{self, AsyncWrite};

pub use self::builder::Builder;
use self::{deflate::Deflate, deflater::Deflater};
use crate::writer::MAX_BUF_SIZE;

#[cfg(feature = "libdeflate")]
type CompressionLevel = libdeflater::CompressionLvl;
#[cfg(not(feature = "libdeflate"))]
type CompressionLevel = flate2::Compression;

pin_project! {
    /// An async BGZF writer.
    pub struct Writer<W> {
        #[pin]
        sink: Buffer<Deflater<W>, Deflate>,
        buf: BytesMut,
        #[pin]
        eof_buf: Bytes,
        compression_level: CompressionLevel,
    }
}

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Creates an async BGZF writer with a default compression level.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::r#async::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Builder::default().build_from_writer(inner)
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let writer = bgzf::r#async::Writer::new(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> W {
        self.sink.into_inner().into_inner()
    }

    fn remaining(&self) -> usize {
        MAX_BUF_SIZE - self.buf.len()
    }

    fn has_remaining(&self) -> bool {
        self.buf.len() < MAX_BUF_SIZE
    }
}

impl<W> AsyncWrite for Writer<W>
where
    W: AsyncWrite + Unpin,
{
    fn poll_write(
        mut self: Pin<&mut Self>,
        cx: &mut Context<'_>,
        buf: &[u8],
    ) -> Poll<io::Result<usize>> {
        if !self.has_remaining() {
            if let Err(e) = ready!(self.as_mut().poll_flush(cx)) {
                return Poll::Ready(Err(e));
            }
        }

        let amt = self.remaining().min(buf.len());
        self.as_mut().buf.extend_from_slice(&buf[..amt]);

        Poll::Ready(Ok(amt))
    }

    fn poll_flush(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<io::Result<()>> {
        let mut this = self.project();

        if this.buf.is_empty() {
            return Poll::Ready(Ok(()));
        }

        ready!(this.sink.as_mut().poll_ready(cx))?;

        let buf = this.buf.split();
        this.sink
            .as_mut()
            .start_send(Deflate::new(buf, *this.compression_level))?;

        Poll::Ready(Ok(()))
    }

    fn poll_shutdown(mut self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<io::Result<()>> {
        ready!(self.as_mut().poll_flush(cx))?;

        let mut this = self.project();
        let mut sink = this.sink.as_mut();

        ready!(sink.as_mut().poll_close(cx))?;

        let mut inner = sink.get_mut().get_mut().get_mut();

        while this.eof_buf.has_remaining() {
            let bytes_written = ready!(Pin::new(&mut inner).poll_write(cx, this.eof_buf.chunk()))?;

            this.eof_buf.advance(bytes_written);

            if bytes_written == 0 {
                return Poll::Ready(Err(io::Error::from(io::ErrorKind::WriteZero)));
            }
        }

        Poll::Ready(Ok(()))
    }
}
