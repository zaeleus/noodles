use std::{
    pin::Pin,
    task::{Context, Poll, ready},
    vec,
};

use noodles_bgzf as bgzf;
use pin_project_lite::pin_project;
use tokio::io::{self, AsyncBufRead, AsyncRead, AsyncSeek, ReadBuf};

use crate::binning_index::index::reference_sequence::bin::Chunk;

enum State {
    Seek,
    Seeking(Chunk),
    Read(bgzf::VirtualPosition),
    Done,
}

pin_project! {
    /// An async query reader.
    ///
    /// This reader returns the uncompressed data between all the given chunks.
    pub struct Query<'r, R>
    where
        R: AsyncRead,
    {
        #[pin]
        reader: &'r mut bgzf::r#async::io::Reader<R>,
        chunks: vec::IntoIter<Chunk>,
        state: State,
    }
}

impl<'r, R> Query<'r, R>
where
    R: AsyncRead + Unpin,
{
    /// Creates an async query reader.
    pub fn new(reader: &'r mut bgzf::r#async::io::Reader<R>, chunks: Vec<Chunk>) -> Self {
        Self {
            reader,
            chunks: chunks.into_iter(),
            state: State::Seek,
        }
    }
}

impl<R> AsyncRead for Query<'_, R>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    fn poll_read(
        mut self: Pin<&mut Self>,
        cx: &mut Context<'_>,
        buf: &mut ReadBuf<'_>,
    ) -> Poll<io::Result<()>> {
        let src = ready!(self.as_mut().poll_fill_buf(cx))?;

        let amt = src.len().min(buf.remaining());
        buf.put_slice(&src[..amt]);

        self.consume(amt);

        Poll::Ready(Ok(()))
    }
}

impl<R> AsyncBufRead for Query<'_, R>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    fn poll_fill_buf(self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<io::Result<&[u8]>> {
        let mut this = self.project();

        loop {
            *this.state = match this.state {
                State::Seek => match this.chunks.next() {
                    Some(chunk) => State::Seeking(chunk),
                    None => State::Done,
                },
                State::Seeking(chunk) => {
                    ready!(this.reader.as_mut().poll_seek(cx, chunk.start()))?;
                    State::Read(chunk.end())
                }
                State::Read(chunk_end) => {
                    if this.reader.as_ref().virtual_position() < *chunk_end {
                        return this.reader.poll_fill_buf(cx);
                    } else {
                        State::Seek
                    }
                }
                State::Done => return Poll::Ready(Ok(&[])),
            };
        }
    }

    fn consume(self: Pin<&mut Self>, amt: usize) {
        self.project().reader.consume(amt);
    }
}
