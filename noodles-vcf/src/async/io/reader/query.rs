use std::vec;

use futures::{Stream, stream};
use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::Reader;
use crate::{Header, Record, io::reader::query::intersects};

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    Done,
}

struct Query<'r, 'h: 'r, R>
where
    R: AsyncRead + AsyncSeek,
{
    inner: &'r mut Reader<bgzf::r#async::io::Reader<R>>,
    chunks: vec::IntoIter<Chunk>,

    header: &'h Header,
    reference_sequence_name: Vec<u8>,
    interval: Interval,

    state: State,
}

impl<'r, 'h: 'r, R> Query<'r, 'h, R>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    fn new(
        inner: &'r mut Reader<bgzf::r#async::io::Reader<R>>,
        chunks: Vec<Chunk>,
        header: &'h Header,
        reference_sequence_name: Vec<u8>,
        interval: Interval,
    ) -> Self {
        Self {
            inner,
            chunks: chunks.into_iter(),
            header,
            reference_sequence_name,
            interval,
            state: State::Seek,
        }
    }

    async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        loop {
            match self.state {
                State::Seek => {
                    self.state = match self.chunks.next() {
                        Some(chunk) => {
                            self.inner.get_mut().seek(chunk.start()).await?;
                            State::Read(chunk.end())
                        }
                        None => State::Done,
                    };
                }
                State::Read(chunk_end) => match self.inner.read_record(record).await? {
                    0 => self.state = State::Seek,
                    n => {
                        if self.inner.get_ref().virtual_position() >= chunk_end {
                            self.state = State::Seek;
                        }

                        if intersects(
                            self.header,
                            record,
                            &self.reference_sequence_name,
                            self.interval,
                        )? {
                            return Ok(n);
                        }
                    }
                },
                State::Done => return Ok(0),
            }
        }
    }
}

pub fn query<'r, 'h: 'r, R>(
    inner: &'r mut Reader<bgzf::r#async::io::Reader<R>>,
    chunks: Vec<Chunk>,
    header: &'h Header,
    reference_sequence_name: Vec<u8>,
    interval: Interval,
) -> impl Stream<Item = io::Result<Record>> + 'r
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    let ctx = Query::new(inner, chunks, header, reference_sequence_name, interval);

    Box::pin(stream::try_unfold(ctx, |mut ctx| async {
        let mut record = Record::default();

        match ctx.read_record(&mut record).await? {
            0 => Ok(None),
            _ => Ok(Some((record, ctx))),
        }
    }))
}
