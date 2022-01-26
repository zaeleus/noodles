use std::ops::RangeBounds;

use futures::{stream, Stream};
use noodles_bgzf as bgzf;
use noodles_csi::index::reference_sequence::bin::Chunk;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::Reader;
use crate::{
    reader::query::{intersects, next_chunk, resolve_interval},
    Header, Record,
};

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    Done,
}

struct Context<'r, R>
where
    R: AsyncRead + AsyncSeek,
{
    reader: &'r mut Reader<bgzf::AsyncReader<R>>,

    chunks: Vec<Chunk>,
    i: usize,

    reference_sequence_name: String,
    start: i32,
    end: i32,

    state: State,

    header: &'r Header,
}

pub fn query<'r, R, B>(
    reader: &'r mut Reader<bgzf::AsyncReader<R>>,
    chunks: Vec<Chunk>,
    reference_sequence_name: String,
    interval: B,
    header: &'r Header,
) -> impl Stream<Item = io::Result<Record>> + 'r
where
    R: AsyncRead + AsyncSeek + Unpin,
    B: RangeBounds<i32>,
{
    let (start, end) = resolve_interval(interval);

    let ctx = Context {
        reader,

        chunks,
        i: 0,

        reference_sequence_name,
        start,
        end,

        state: State::Seek,

        header,
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
                State::Read(chunk_end) => match next_record(ctx.reader, ctx.header).await? {
                    Some(record) => {
                        if ctx.reader.virtual_position() >= chunk_end {
                            ctx.state = State::Seek;
                        }

                        if intersects(&record, &ctx.reference_sequence_name, ctx.start, ctx.end)? {
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

async fn next_record<R>(
    reader: &mut Reader<bgzf::AsyncReader<R>>,
    header: &Header,
) -> io::Result<Option<Record>>
where
    R: AsyncRead + Unpin,
{
    let mut buf = String::new();

    match reader.read_record(&mut buf).await? {
        0 => Ok(None),
        _ => Record::try_from_str(&buf, header)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}
