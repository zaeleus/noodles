//! Async BGZF reader.

mod builder;
mod inflate;
mod inflater;

use std::{
    cmp, io,
    num::NonZero,
    pin::Pin,
    task::{Context, Poll, ready},
};

use futures::{Stream, TryStreamExt, stream::TryBuffered};
use pin_project_lite::pin_project;
use tokio::io::{AsyncBufRead, AsyncRead, AsyncSeek, ReadBuf};

pub use self::builder::Builder;
use self::inflater::Inflater;
use crate::{
    VirtualPosition, gzi,
    io::{Block, reader::frame::block_initialize},
};

enum SeekState<R>
where
    R: AsyncRead,
{
    Init,
    Seek(TryBuffered<Inflater<R>>),
    Finish(TryBuffered<Inflater<R>>),
    Done,
}

pin_project! {
    /// An async BGZF reader.
    pub struct Reader<R>
    where
        R: AsyncRead,
    {
        #[pin]
        stream: Option<TryBuffered<Inflater<R>>>,
        block: Block,
        position: u64,
        worker_count: NonZero<usize>,
        seek_state: Option<SeekState<R>>,
    }
}

impl<R> Reader<R>
where
    R: AsyncRead,
{
    /// Creates an async BGZF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use tokio::io;
    /// let reader = bgzf::r#async::io::Reader::new(io::empty());
    /// ```
    pub fn new(inner: R) -> Self {
        Builder::default().build_from_reader(inner)
    }

    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use tokio::io;
    /// let reader = bgzf::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &R {
        let stream = self.stream.as_ref().expect("missing stream");
        stream.get_ref().get_ref()
    }

    /// Returns a mutable reference to the underlying stream.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use tokio::io;
    /// let mut reader = bgzf::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        let stream = self.stream.as_mut().expect("missing stream");
        stream.get_mut().get_mut()
    }

    /// Returns a pinned mutable reference to the underlying stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::pin::Pin;
    /// use noodles_bgzf as bgzf;
    /// use tokio::io;
    /// let mut reader = bgzf::r#async::io::Reader::new(io::empty());
    /// let _inner = Pin::new(&mut reader).get_pin_mut();
    /// ```
    pub fn get_pin_mut(self: Pin<&mut Self>) -> Pin<&mut R> {
        let stream = self.project().stream.as_pin_mut().expect("missing stream");
        stream.get_pin_mut().get_pin_mut()
    }

    /// Unwraps and returns the underlying stream.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use tokio::io;
    /// let reader = bgzf::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> R {
        let stream = self.stream.expect("missing stream");
        stream.into_inner().into_inner()
    }

    /// Returns the current position of the stream.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use tokio::io;
    /// let reader = bgzf::r#async::io::Reader::new(io::empty());
    /// assert_eq!(reader.position(), 0);
    /// ```
    pub fn position(&self) -> u64 {
        self.position
    }

    /// Returns the current virtual position of the stream.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// use tokio::io;
    /// let reader = bgzf::r#async::io::Reader::new(io::empty());
    /// assert_eq!(reader.virtual_position(), bgzf::VirtualPosition::from(0));
    /// ```
    pub fn virtual_position(&self) -> VirtualPosition {
        self.block.virtual_position()
    }
}

impl<R> Reader<R>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    /// Seeks the stream to the given virtual position.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_bgzf as bgzf;
    /// use tokio::io;
    /// let mut reader = bgzf::r#async::io::Reader::new(io::empty());
    /// reader.seek(bgzf::VirtualPosition::MIN).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn seek(&mut self, pos: VirtualPosition) -> io::Result<VirtualPosition> {
        let (cpos, upos) = pos.into();

        let stream = self.stream.as_mut().expect("missing stream");
        let blocks = stream.get_mut();

        blocks.seek(pos).await?;

        let stream = self.stream.take().expect("missing stream");
        let blocks = stream.into_inner();

        let mut stream = blocks.try_buffered(self.worker_count.get());
        let item = stream.try_next().await;
        self.stream.replace(stream);

        match item {
            Ok(Some(block)) => {
                self.position = cpos + block.size();
                self.block = block;
                self.block.set_position(cpos);
            }
            Ok(None) => {
                self.position = cpos;
                block_initialize(&mut self.block, 0, 0);
                self.block.set_position(cpos);
            }
            Err(e) => {
                self.position = cpos;
                block_initialize(&mut self.block, 0, 0);
                self.block.set_position(cpos);
                return Err(e);
            }
        };

        let data = self.block.data_mut();

        if usize::from(upos) <= data.len() {
            data.set_position(usize::from(upos));
        } else {
            return Err(io::Error::from(io::ErrorKind::InvalidInput));
        }

        Ok(pos)
    }

    #[doc(hidden)]
    pub fn poll_seek(
        mut self: Pin<&mut &mut Self>,
        cx: &mut Context<'_>,
        pos: VirtualPosition,
    ) -> Poll<io::Result<VirtualPosition>> {
        loop {
            self.seek_state = match self.seek_state.take().unwrap() {
                SeekState::Init => {
                    let stream = self.stream.take().expect("missing stream");
                    Some(SeekState::Seek(stream))
                }
                SeekState::Seek(mut stream) => {
                    let blocks = Pin::new(&mut stream).get_pin_mut();

                    match blocks.poll_seek(cx, pos) {
                        Poll::Ready(Ok(_)) => {
                            let blocks = stream.into_inner();
                            let stream = blocks.try_buffered(self.worker_count.get());
                            Some(SeekState::Finish(stream))
                        }
                        Poll::Ready(Err(e)) => {
                            self.stream.replace(stream);
                            self.seek_state = Some(SeekState::Init);
                            return Poll::Ready(Err(e));
                        }
                        Poll::Pending => {
                            self.seek_state = Some(SeekState::Seek(stream));
                            return Poll::Pending;
                        }
                    }
                }
                SeekState::Finish(mut stream) => {
                    let (cpos, upos) = pos.into();

                    let item = match Pin::new(&mut stream).poll_next(cx) {
                        Poll::Ready(item) => item,
                        Poll::Pending => {
                            self.seek_state = Some(SeekState::Finish(stream));
                            return Poll::Pending;
                        }
                    };

                    self.stream.replace(stream);

                    match item {
                        Some(Ok(block)) => {
                            self.position = cpos + block.size();
                            self.block = block;
                            self.block.set_position(cpos);
                        }
                        Some(Err(e)) => {
                            self.position = cpos;
                            block_initialize(&mut self.block, 0, 0);
                            self.block.set_position(cpos);

                            self.seek_state = Some(SeekState::Init);

                            return Poll::Ready(Err(e));
                        }
                        None => {
                            self.position = cpos;
                            block_initialize(&mut self.block, 0, 0);
                            self.block.set_position(cpos);
                        }
                    };

                    let data = self.block.data_mut();

                    if usize::from(upos) <= data.len() {
                        data.set_position(usize::from(upos));
                    } else {
                        self.seek_state = Some(SeekState::Init);
                        return Poll::Ready(Err(io::Error::from(io::ErrorKind::InvalidInput)));
                    }

                    self.seek_state = Some(SeekState::Done);

                    return Poll::Ready(Ok(pos));
                }
                SeekState::Done => Some(SeekState::Init),
            };
        }
    }

    /// Seeks the stream to the given uncompressed position.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_bgzf::{self as bgzf, gzi};
    /// use tokio::io;
    ///
    /// let mut reader = bgzf::r#async::io::Reader::new(io::empty());
    ///
    /// let index = gzi::Index::default();
    /// reader.seek_by_uncompressed_position(&index, 0).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn seek_by_uncompressed_position(
        &mut self,
        index: &gzi::Index,
        pos: u64,
    ) -> io::Result<u64> {
        let virtual_position = index.query(pos)?;
        self.seek(virtual_position).await?;
        Ok(pos)
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
        let src = ready!(self.as_mut().poll_fill_buf(cx))?;

        let amt = cmp::min(src.len(), buf.remaining());
        buf.put_slice(&src[..amt]);

        self.consume(amt);

        Poll::Ready(Ok(()))
    }
}

impl<R> AsyncBufRead for Reader<R>
where
    R: AsyncRead,
{
    fn poll_fill_buf(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<io::Result<&[u8]>> {
        let this = self.project();

        if !this.block.data().has_remaining() {
            let mut stream = this.stream.as_pin_mut().expect("missing stream");

            loop {
                match ready!(stream.as_mut().poll_next(cx)) {
                    Some(Ok(mut block)) => {
                        block.set_position(*this.position);
                        *this.position += block.size();
                        let data_len = block.data().len();
                        *this.block = block;

                        if data_len > 0 {
                            break;
                        }
                    }
                    Some(Err(e)) => return Poll::Ready(Err(e)),
                    None => {
                        block_initialize(this.block, 0, 0);
                        this.block.set_position(*this.position);
                        break;
                    }
                }
            }
        }

        Poll::Ready(Ok(this.block.data().as_ref()))
    }

    fn consume(self: Pin<&mut Self>, amt: usize) {
        let this = self.project();
        this.block.data_mut().consume(amt);
    }
}

#[cfg(test)]
mod tests {
    use std::{future, io::Cursor};

    use tokio::io::AsyncReadExt;

    use super::*;

    #[tokio::test]
    async fn test_read_with_empty_block() -> io::Result<()> {
        #[rustfmt::skip]
        let data = [
            // block 0 (b"noodles")
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x22, 0x00, 0xcb, 0xcb, 0xcf, 0x4f, 0xc9, 0x49, 0x2d, 0x06, 0x00, 0xa1,
            0x58, 0x2a, 0x80, 0x07, 0x00, 0x00, 0x00,
            // block 1 (b"")
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            // block 2 (b"bgzf")
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1f, 0x00, 0x4b, 0x4a, 0xaf, 0x4a, 0x03, 0x00, 0x20, 0x68, 0xf2, 0x8c,
            0x04, 0x00, 0x00, 0x00,
            // EOF block
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];

        let mut reader = Reader::new(&data[..]);
        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).await?;

        assert_eq!(buf, b"noodlesbgzf");

        Ok(())
    }

    #[tokio::test]
    async fn test_seek() -> Result<(), Box<dyn std::error::Error>> {
        #[rustfmt::skip]
        let data = [
            // block 0, udata = b"noodles"
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x22, 0x00, 0xcb, 0xcb, 0xcf, 0x4f, 0xc9, 0x49, 0x2d, 0x06, 0x00, 0xa1,
            0x58, 0x2a, 0x80, 0x07, 0x00, 0x00, 0x00,
            // EOF block
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];

        let mut reader = Reader::new(Cursor::new(&data));

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).await?;

        let eof = VirtualPosition::try_from((63, 0))?;
        assert_eq!(reader.virtual_position(), eof);

        let position = VirtualPosition::try_from((0, 3))?;
        reader.seek(position).await?;

        assert_eq!(reader.virtual_position(), position);

        buf.clear();
        reader.read_to_end(&mut buf).await?;

        assert_eq!(buf, b"dles");
        assert_eq!(reader.virtual_position(), eof);

        Ok(())
    }

    #[tokio::test]
    async fn test_seek_with_uncompressed_position_gt_data_len()
    -> Result<(), crate::virtual_position::TryFromU64U16TupleError> {
        #[rustfmt::skip]
        let data = [
            // block 0 (b"noodles")
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x22, 0x00, 0xcb, 0xcb, 0xcf, 0x4f, 0xc9, 0x49, 0x2d, 0x06, 0x00, 0xa1,
            0x58, 0x2a, 0x80, 0x07, 0x00, 0x00, 0x00,
            // EOF block
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];

        let mut reader = Reader::new(Cursor::new(&data));

        assert!(matches!(
            reader.seek(VirtualPosition::try_from((0, 8))?).await,
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[tokio::test]
    async fn test_poll_seek() -> io::Result<()> {
        #[rustfmt::skip]
        let data = [
            // block 0 (b"noodles")
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x22, 0x00, 0xcb, 0xcb, 0xcf, 0x4f, 0xc9, 0x49, 0x2d, 0x06, 0x00, 0xa1,
            0x58, 0x2a, 0x80, 0x07, 0x00, 0x00, 0x00,
            // EOF block
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];

        let mut reader = Reader::new(Cursor::new(&data));

        future::poll_fn(|cx| Pin::new(&mut &mut reader).poll_seek(cx, VirtualPosition::MIN))
            .await?;

        assert_eq!(reader.virtual_position(), VirtualPosition::MIN);

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).await?;
        assert_eq!(buf, b"noodles");

        future::poll_fn(|cx| Pin::new(&mut &mut reader).poll_seek(cx, VirtualPosition::MIN))
            .await?;

        assert_eq!(reader.virtual_position(), VirtualPosition::MIN);

        buf.clear();
        reader.read_to_end(&mut buf).await?;
        assert_eq!(buf, b"noodles");

        Ok(())
    }
}
