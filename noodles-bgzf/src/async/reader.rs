mod block_decoder;
mod blocks;

use std::{
    cmp, io,
    pin::Pin,
    task::{Context, Poll},
};

use futures::{stream::TryBuffered, Stream, TryStreamExt};
use pin_project_lite::pin_project;
use tokio::io::{AsyncBufRead, AsyncRead, ReadBuf};

use crate::{Block, VirtualPosition};

use self::blocks::Blocks;

const WORKER_COUNT: usize = 8;

pin_project! {
    /// An async BGZF reader.
    pub struct Reader<R>
    where
        R: AsyncRead,
    {
        #[pin]
        stream: TryBuffered<Blocks<R>>,
        block: Block,
    }
}

impl<R> Reader<R>
where
    R: AsyncRead,
{
    /// Creates an async BGZF reader.
    pub fn new(inner: R) -> Self {
        Self {
            stream: Blocks::new(inner).try_buffered(WORKER_COUNT),
            block: Block::default(),
        }
    }

    /// Returns the current virtual position of the stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// let data = [];
    /// let reader = bgzf::AsyncReader::new(&data[..]);
    /// assert_eq!(reader.virtual_position(), bgzf::VirtualPosition::from(0));
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn virtual_position(&self) -> VirtualPosition {
        self.block.virtual_position()
    }
}

impl<R> AsyncRead for Reader<R>
where
    R: AsyncRead,
{
    fn poll_read(
        mut self: Pin<&mut Self>,
        cx: &mut Context<'_>,
        buf: &mut ReadBuf<'_>,
    ) -> Poll<io::Result<()>> {
        let block_buf = match self.as_mut().poll_fill_buf(cx) {
            Poll::Ready(Ok(b)) => b,
            Poll::Ready(Err(e)) => return Poll::Ready(Err(e)),
            Poll::Pending => return Poll::Pending,
        };

        let len = cmp::min(block_buf.len(), buf.remaining());
        buf.put_slice(&block_buf[..len]);

        self.consume(len);

        Poll::Ready(Ok(()))
    }
}

impl<R> AsyncBufRead for Reader<R>
where
    R: AsyncRead,
{
    fn poll_fill_buf(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<io::Result<&[u8]>> {
        let this = self.project();

        if this.block.is_eof() {
            match this.stream.poll_next(cx) {
                Poll::Ready(Some(Ok(block))) => {
                    *this.block = block;
                }
                Poll::Ready(Some(Err(e))) => return Poll::Ready(Err(e)),
                Poll::Ready(None) => return Poll::Ready(Ok(&[])),
                Poll::Pending => return Poll::Pending,
            }
        }

        return Poll::Ready(Ok(this.block.fill_buf()));
    }

    fn consume(self: Pin<&mut Self>, amt: usize) {
        let this = self.project();
        let upos = cmp::min(this.block.ulen(), this.block.upos() + amt as u32);
        this.block.set_upos(upos);
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use tokio::io::AsyncReadExt;

    use crate::Writer;

    use super::*;

    #[tokio::test]
    async fn test_read() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());
        writer.write_all(b"noodles")?;

        let data = writer.finish()?;
        let mut reader = Reader::new(&data[..]);

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).await?;

        assert_eq!(buf, b"noodles");

        Ok(())
    }
}
