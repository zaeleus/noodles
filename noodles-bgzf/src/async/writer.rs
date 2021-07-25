pub(crate) mod deflate;
mod deflater;

use std::{
    cmp,
    future::Future,
    pin::Pin,
    task::{Context, Poll},
};

use bytes::BytesMut;
use futures::{ready, sink::Buffer, Sink, SinkExt};
use pin_project_lite::pin_project;
use tokio::{
    io::{self, AsyncWrite},
    pin,
};
use tokio_util::codec::FramedWrite;

use super::BlockCodec;

use crate::{block, writer::BGZF_EOF};

use self::{deflate::Deflate, deflater::Deflater};

const WORKER_COUNT: usize = 8;

pin_project! {
    /// An async BGZF writer.
    pub struct Writer<W> {
        #[pin]
        sink: Buffer<Deflater<W>, Deflate>,
        buf: BytesMut,
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
    /// let writer = bgzf::AsyncWriter::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self {
            sink: Deflater::new(FramedWrite::new(inner, BlockCodec)).buffer(WORKER_COUNT),
            buf: BytesMut::with_capacity(block::MAX_UNCOMPRESSED_DATA_LENGTH),
        }
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
    ) -> Poll<Result<usize, io::Error>> {
        if self.buf.len() >= block::MAX_UNCOMPRESSED_DATA_LENGTH {
            if let Err(e) = ready!(self.as_mut().poll_flush(cx)) {
                return Poll::Ready(Err(e));
            }
        }

        let n = cmp::min(
            block::MAX_UNCOMPRESSED_DATA_LENGTH - self.buf.len(),
            buf.len(),
        );

        self.as_mut().buf.extend_from_slice(&buf[..n]);

        Poll::Ready(Ok(n))
    }

    fn poll_flush(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Result<(), io::Error>> {
        let mut this = self.project();

        if this.buf.is_empty() {
            return Poll::Ready(Ok(()));
        }

        ready!(this.sink.as_mut().poll_ready(cx))?;

        let buf = this.buf.split();
        this.sink.as_mut().start_send(Deflate::new(buf))?;

        Poll::Ready(Ok(()))
    }

    fn poll_shutdown(
        mut self: Pin<&mut Self>,
        cx: &mut Context<'_>,
    ) -> Poll<Result<(), io::Error>> {
        use tokio::io::AsyncWriteExt;

        let mut this = self.as_mut();

        ready!(this.as_mut().poll_flush(cx))?;

        let mut sink = this.project().sink;
        ready!(sink.as_mut().poll_close(cx))?;

        let write_all = sink.get_mut().get_mut().get_mut().write_all(BGZF_EOF);
        pin!(write_all);
        write_all.poll(cx)
    }
}
