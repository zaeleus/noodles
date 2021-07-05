use std::{
    future::Future,
    io::{self, Read},
    pin::Pin,
    task::{Context, Poll},
};

use bytes::BytesMut;
use flate2::bufread::DeflateDecoder;
use futures::Stream;
use pin_project_lite::pin_project;
use tokio::io::AsyncRead;
use tokio_util::codec::FramedRead;

use crate::{gz, Block, BGZF_HEADER_SIZE};

use super::block_decoder::BlockDecoder;

pin_project! {
    pub struct Blocks<R> {
        #[pin]
        inner: FramedRead<R, BlockDecoder>,
    }
}

impl<R> Blocks<R>
where
    R: AsyncRead,
{
    pub fn new(inner: R) -> Self {
        Self {
            inner: FramedRead::new(inner, BlockDecoder),
        }
    }
}

impl<R> Stream for Blocks<R>
where
    R: AsyncRead,
{
    type Item = io::Result<Pin<Box<dyn Future<Output = io::Result<Block>>>>>;

    fn poll_next(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        match self.project().inner.poll_next(cx) {
            Poll::Ready(Some(Ok(buf))) => Poll::Ready(Some(Ok(Box::pin(build_block(buf))))),
            Poll::Ready(Some(Err(e))) => Poll::Ready(Some(Err(e))),
            Poll::Ready(None) => Poll::Ready(None),
            Poll::Pending => Poll::Pending,
        }
    }
}

async fn build_block(src: BytesMut) -> io::Result<Block> {
    tokio::task::spawn_blocking(move || {
        let cdata_len = src.len() - BGZF_HEADER_SIZE - gz::TRAILER_SIZE;
        let cdata = &src[BGZF_HEADER_SIZE..BGZF_HEADER_SIZE + cdata_len];

        let mut block = Block::default();
        let udata = block.data_mut();

        let mut decoder = DeflateDecoder::new(cdata);
        decoder.read_to_end(udata)?;

        Ok(block)
    })
    .await?
}
