use std::{
    future::Future,
    io::{self, Read},
    mem,
    pin::Pin,
    task::{Context, Poll},
};

use bytes::{Buf, BytesMut};
use flate2::bufread::DeflateDecoder;
use futures::{ready, Stream};
use pin_project_lite::pin_project;
use tokio::io::{AsyncRead, AsyncSeek, AsyncSeekExt, SeekFrom};
use tokio_util::codec::FramedRead;

use crate::{gz, Block, VirtualPosition, BGZF_HEADER_SIZE};

use crate::r#async::BlockCodec;

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
    type Item = io::Result<Pin<Box<dyn Future<Output = io::Result<Block>> + Send>>>;

    fn poll_next(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        match ready!(self.project().inner.poll_next(cx)) {
            Some(Ok(buf)) => Poll::Ready(Some(Ok(Box::pin(inflate(buf))))),
            Some(Err(e)) => Poll::Ready(Some(Err(e))),
            None => Poll::Ready(None),
        }
    }
}

async fn inflate(mut src: BytesMut) -> io::Result<Block> {
    tokio::task::spawn_blocking(move || {
        let mut header = src.split_to(BGZF_HEADER_SIZE);
        header.advance(16); // [ID1, ..., SLEN]
        let bsize = u64::from(header.get_u16_le()) + 1;

        let cdata = src.split_to(src.len() - gz::TRAILER_SIZE);

        // trailer
        src.advance(mem::size_of::<u32>()); // CRC32
        let r#isize = src.get_u32_le();

        let mut block = Block::default();
        let udata = block.data_mut();
        udata.reserve_exact((r#isize + 1) as usize);

        let mut decoder = DeflateDecoder::new(&cdata[..]);
        decoder.read_to_end(udata)?;

        block.set_clen(bsize);

        Ok(block)
    })
    .await?
}
