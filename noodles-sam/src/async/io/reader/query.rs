use std::vec;

use futures::{stream, Stream};
use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::{Header, Reader};
use crate::{io::reader::query::intersects, Record};

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    Done,
}

struct Context<'r, 'h: 'r, R>
where
    R: AsyncRead + AsyncSeek,
{
    reader: &'r mut Reader<bgzf::AsyncReader<R>>,
    chunks: vec::IntoIter<Chunk>,

    header: &'h Header,
    reference_sequence_id: usize,
    interval: Interval,

    state: State,
}

pub(super) fn query<'r, 'h: 'r, R>(
    reader: &'r mut Reader<bgzf::AsyncReader<R>>,
    chunks: Vec<Chunk>,
    header: &'h Header,
    reference_sequence_id: usize,
    interval: Interval,
) -> impl Stream<Item = io::Result<Record>> + 'r
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    let ctx = Context {
        reader,
        chunks: chunks.into_iter(),

        header,
        reference_sequence_id,
        interval,

        state: State::Seek,
    };

    Box::pin(stream::try_unfold(ctx, |mut ctx| async {
        loop {
            match ctx.state {
                State::Seek => {
                    ctx.state = match ctx.chunks.next() {
                        Some(chunk) => {
                            ctx.reader.get_mut().seek(chunk.start()).await?;
                            State::Read(chunk.end())
                        }
                        None => State::Done,
                    };
                }
                State::Read(chunk_end) => match next_record(ctx.reader).await? {
                    Some(record) => {
                        if ctx.reader.get_ref().virtual_position() >= chunk_end {
                            ctx.state = State::Seek;
                        }

                        if intersects(ctx.header, &record, ctx.reference_sequence_id, ctx.interval)?
                        {
                            return Ok(Some((record, ctx)));
                        }
                    }
                    None => ctx.state = State::Seek,
                },
                State::Done => return Ok(None),
            }
        }
    }))
}

async fn next_record<R>(reader: &mut Reader<bgzf::AsyncReader<R>>) -> io::Result<Option<Record>>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    let mut record = Record::default();

    reader.read_record(&mut record).await.map(|n| match n {
        0 => None,
        _ => Some(record),
    })
}
