use std::{
    pin::Pin,
    task::{Context, Poll},
};

use flate2::Crc;
use futures::ready;
use pin_project_lite::pin_project;
use tokio::io::{self, AsyncRead, ReadBuf};

pin_project! {
    #[derive(Debug)]
    pub struct CrcReader<R> {
        #[pin]
        inner: R,
        crc: Crc,
    }
}

impl<R> CrcReader<R>
where
    R: AsyncRead + Unpin,
{
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            crc: Crc::new(),
        }
    }

    pub fn crc(&self) -> &Crc {
        &self.crc
    }

    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> AsyncRead for CrcReader<R>
where
    R: AsyncRead + Unpin,
{
    fn poll_read(
        self: Pin<&mut Self>,
        cx: &mut Context<'_>,
        buf: &mut ReadBuf<'_>,
    ) -> Poll<io::Result<()>> {
        let this = self.project();
        ready!(this.inner.poll_read(cx, buf))?;
        this.crc.update(buf.filled());
        Poll::Ready(Ok(()))
    }
}
