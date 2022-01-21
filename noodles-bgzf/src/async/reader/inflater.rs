use std::{
    io,
    pin::Pin,
    task::{Context, Poll},
};

use futures::{ready, Stream};
use pin_project_lite::pin_project;
use tokio::io::{AsyncRead, AsyncSeek, AsyncSeekExt, SeekFrom};
use tokio_util::codec::FramedRead;

use super::inflate::Inflate;
use crate::{r#async::BlockCodec, VirtualPosition};

pin_project! {
    pub struct Inflater<R> {
        #[pin]
        inner: FramedRead<R, BlockCodec>,
    }
}

impl<R> Inflater<R>
where
    R: AsyncRead,
{
    pub fn new(inner: R) -> Self {
        Self {
            inner: FramedRead::new(inner, BlockCodec),
        }
    }
}

impl<R> Inflater<R>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    pub async fn seek(&mut self, pos: VirtualPosition) -> io::Result<VirtualPosition> {
        let cpos = pos.compressed();
        self.inner.get_mut().seek(SeekFrom::Start(cpos)).await?;

        self.inner.read_buffer_mut().clear();

        Ok(pos)
    }
}

impl<R> Stream for Inflater<R>
where
    R: AsyncRead,
{
    type Item = io::Result<Inflate>;

    fn poll_next(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        match ready!(self.project().inner.poll_next(cx)) {
            Some(Ok(buf)) => Poll::Ready(Some(Ok(Inflate::new(buf)))),
            Some(Err(e)) => Poll::Ready(Some(Err(e))),
            None => Poll::Ready(None),
        }
    }
}
