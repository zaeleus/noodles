use std::ops::RangeBounds;

use futures::{stream, Stream};
use noodles_bgzf as bgzf;
use noodles_csi::index::reference_sequence::bin::Chunk;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::Reader;
use crate::{
    reader::query::{intersects, next_chunk, resolve_interval},
    Record,
};

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    Done,
}

struct Context<'a, R>
where
    R: AsyncRead + AsyncSeek,
{
    reader: &'a mut Reader<R>,

    chunks: Vec<Chunk>,
    i: usize,

    reference_sequence_id: usize,
    start: i32,
    end: i32,

    state: State,
}

pub fn query<R, B>(
    reader: &mut Reader<R>,
    chunks: Vec<Chunk>,
    reference_sequence_id: usize,
    interval: B,
) -> impl Stream<Item = io::Result<Record>> + '_
where
    R: AsyncRead + AsyncSeek + Unpin,
    B: RangeBounds<i32>,
{
    let (start, end) = resolve_interval(interval);

    let ctx = Context {
        reader,

        chunks,
        i: 0,

        reference_sequence_id,
        start,
        end,

        state: State::Seek,
    };

    Box::pin(stream::try_unfold(ctx, |mut ctx| async {
        loop {
            match ctx.state {
                State::Seek => {
                    ctx.state = match next_chunk(&ctx.chunks, &mut ctx.i) {
                        Some(chunk) => {
                            ctx.reader.seek(chunk.start()).await?;
                            State::Read(chunk.end())
                        }
                        None => State::Done,
                    };
                }
                State::Read(chunk_end) => match next_record(&mut ctx.reader).await? {
                    Some(record) => {
                        if ctx.reader.virtual_position() >= chunk_end {
                            ctx.state = State::Seek;
                        }

                        if intersects(&record, ctx.reference_sequence_id, ctx.start, ctx.end)? {
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

async fn next_record<R>(reader: &mut Reader<R>) -> io::Result<Option<Record>>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    let mut record = Record::default();

    reader.read_record(&mut record).await.map(|n| match n {
        0 => None,
        _ => Some(record),
    })
}
