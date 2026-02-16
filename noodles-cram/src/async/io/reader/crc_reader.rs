use std::{
    pin::Pin,
    task::{Context, Poll, ready},
};

use flate2::Crc;
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

    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
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
        let prev_filled = buf.filled().len();
        ready!(this.inner.poll_read(cx, buf))?;
        this.crc.update(&buf.filled()[prev_filled..]);
        Poll::Ready(Ok(()))
    }
}
