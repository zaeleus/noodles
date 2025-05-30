use std::{
    future::Future,
    io,
    pin::Pin,
    task::{Context, Poll, ready},
};

use futures::Sink;
use pin_project_lite::pin_project;
use tokio::io::AsyncWrite;
use tokio_util::codec::FramedWrite;

use super::Deflate;
use crate::r#async::BlockCodec;

pin_project! {
    pub struct Deflater<W> {
        #[pin]
        sink: FramedWrite<W, BlockCodec>,
        #[pin]
        state: Option<Deflate>,
    }
}

impl<W> Deflater<W>
where
    W: AsyncWrite,
{
    pub fn new(sink: FramedWrite<W, BlockCodec>) -> Self {
        Self { sink, state: None }
    }

    pub fn get_mut(&mut self) -> &mut W {
        self.sink.get_mut()
    }

    pub fn into_inner(self) -> W {
        self.sink.into_inner()
    }

    fn poll(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<io::Result<()>> {
        let mut this = self.project();

        let data = match this.state.as_mut().as_pin_mut() {
            Some(deflate) => ready!(deflate.poll(cx))?,
            None => return Poll::Ready(Ok(())),
        };

        this.state.set(None);
        this.sink.start_send(data)?;

        Poll::Ready(Ok(()))
    }
}

impl<W> Sink<Deflate> for Deflater<W>
where
    W: AsyncWrite,
{
    type Error = io::Error;

    fn poll_ready(mut self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Result<(), Self::Error>> {
        ready!(self.as_mut().poll(cx))?;
        self.project().sink.poll_flush(cx)
    }

    fn start_send(self: Pin<&mut Self>, deflate: Deflate) -> Result<(), Self::Error> {
        let mut this = self.project();
        this.state.set(Some(deflate));
        Ok(())
    }

    fn poll_flush(mut self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Result<(), Self::Error>> {
        ready!(self.as_mut().poll(cx))?;
        self.project().sink.poll_flush(cx)
    }

    fn poll_close(mut self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Result<(), Self::Error>> {
        ready!(self.as_mut().poll(cx))?;
        self.project().sink.poll_close(cx)
    }
}
