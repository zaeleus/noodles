use std::ops::RangeBounds;

use futures::{stream, Stream};
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::index::reference_sequence::bin::Chunk;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::Reader;
use crate::{
    reader::query::{intersects, next_chunk},
    Record,
};

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    Done,
}

struct Context<'a, R, B>
where
    R: AsyncRead + AsyncSeek,
    B: RangeBounds<Position> + Copy + 'a,
{
    reader: &'a mut Reader<bgzf::AsyncReader<R>>,

    chunks: Vec<Chunk>,
    i: usize,

    reference_sequence_id: usize,
    interval: B,

    state: State,
}

pub fn query<'a, R, B>(
    reader: &'a mut Reader<bgzf::AsyncReader<R>>,
    chunks: Vec<Chunk>,
    reference_sequence_id: usize,
    interval: B,
) -> impl Stream<Item = io::Result<Record>> + '_
where
    R: AsyncRead + AsyncSeek + Unpin,
    B: RangeBounds<Position> + Copy + 'a,
{
    let ctx = Context {
        reader,

        chunks,
        i: 0,

        reference_sequence_id,
        interval,

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
                State::Read(chunk_end) => match next_record(ctx.reader).await? {
                    Some(record) => {
                        if ctx.reader.virtual_position() >= chunk_end {
                            ctx.state = State::Seek;
                        }

                        if intersects(&record, ctx.reference_sequence_id, ctx.interval)? {
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
