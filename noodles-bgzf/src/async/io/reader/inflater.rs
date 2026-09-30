use std::{
    future, io,
    pin::Pin,
    task::{Context, Poll, ready},
};

use futures::Stream;
use pin_project_lite::pin_project;
use tokio::io::{AsyncRead, AsyncSeek, AsyncSeekExt, ReadBuf, SeekFrom};
use tokio_util::codec::FramedRead;

use super::inflate::Inflate;
use crate::{VirtualPosition, r#async::BlockCodec};

pin_project! {
    pub struct Inflater<R> {
        #[pin]
        inner: FramedRead<R, BlockCodec>,
        is_seeking: bool,
    }
}

impl<R> Inflater<R> {
    pub(super) fn get_ref(&self) -> &R {
        self.inner.get_ref()
    }

    pub(super) fn get_mut(&mut self) -> &mut R {
        self.inner.get_mut()
    }

    pub(super) fn get_pin_mut(self: Pin<&mut Self>) -> Pin<&mut R> {
        self.project().inner.get_pin_mut()
    }

    pub(super) fn into_inner(self) -> R {
        self.inner.into_inner()
    }
}

impl<R> Inflater<R>
where
    R: AsyncRead,
{
    pub fn new(inner: R) -> Self {
        Self {
            inner: FramedRead::new(inner, BlockCodec),
            is_seeking: false,
        }
    }
}

impl<R> Inflater<R>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    pub async fn seek(&mut self, pos: VirtualPosition) -> io::Result<VirtualPosition> {
        // Force read completion.
        let mut buf = ReadBuf::new(&mut []);
        future::poll_fn(|cx| Pin::new(self.inner.get_mut()).poll_read(cx, &mut buf)).await?;

        let cpos = pos.compressed();
        self.inner.get_mut().seek(SeekFrom::Start(cpos)).await?;

        self.inner.read_buffer_mut().clear();

        Ok(pos)
    }

    pub(super) fn poll_seek(
        mut self: Pin<&mut Self>,
        cx: &mut Context<'_>,
        pos: VirtualPosition,
    ) -> Poll<io::Result<VirtualPosition>> {
        let this = self.as_mut().project();
        let mut reader = this.inner.get_pin_mut();

        if !*this.is_seeking {
            // Force read completion.
            let mut buf = ReadBuf::new(&mut []);
            ready!(reader.as_mut().poll_read(cx, &mut buf))?;

            ready!(reader.as_mut().poll_complete(cx))?;

            let cpos = pos.compressed();
            reader.as_mut().start_seek(SeekFrom::Start(cpos))?;

            *this.is_seeking = true;
        }

        match reader.poll_complete(cx) {
            Poll::Ready(Ok(_)) => {}
            Poll::Ready(Err(e)) => {
                *this.is_seeking = false;
                return Poll::Ready(Err(e));
            }
            Poll::Pending => return Poll::Pending,
        }

        *this.is_seeking = false;
        self.inner.read_buffer_mut().clear();

        Poll::Ready(Ok(pos))
    }
}

impl<R> Stream for Inflater<R>
where
    R: AsyncRead,
{
    type Item = io::Result<Inflate>;

    fn poll_next(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        let mut this = self.project();

        match ready!(this.inner.as_mut().poll_next(cx)) {
            Some(Ok(buf)) => Poll::Ready(Some(Ok(Inflate::new(buf)))),
            Some(Err(e)) => Poll::Ready(Some(Err(e))),
            None => match ready!(this.inner.as_mut().poll_next(cx)) {
                Some(Ok(buf)) => Poll::Ready(Some(Ok(Inflate::new(buf)))),
                Some(Err(e)) => Poll::Ready(Some(Err(e))),
                None => Poll::Ready(None),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use futures::StreamExt;

    use super::*;
    use crate::io::writer::BGZF_EOF;

    #[tokio::test]
    async fn test_poll_next() -> io::Result<()> {
        struct R {
            i: usize,
        }

        impl AsyncRead for R {
            fn poll_read(
                mut self: Pin<&mut Self>,
                _cx: &mut Context<'_>,
                buf: &mut ReadBuf<'_>,
            ) -> Poll<io::Result<()>> {
                self.i += 1;

                match self.i {
                    1 => Poll::Ready(Err(io::Error::other(""))),
                    2 => {
                        buf.put_slice(&BGZF_EOF);
                        Poll::Ready(Ok(()))
                    }
                    _ => Poll::Ready(Ok(())),
                }
            }
        }

        let mut inflater = Inflater::new(R { i: 0 });

        assert!(matches!(
            inflater.next().await,
            Some(Err(e)) if e.kind() == io::ErrorKind::Other
        ));

        assert!(matches!(inflater.next().await, Some(Ok(_))));

        assert!(inflater.next().await.is_none());

        Ok(())
    }
}
