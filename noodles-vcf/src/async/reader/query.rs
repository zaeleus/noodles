use std::ops::{Bound, RangeBounds};

use futures::{stream, Stream};
use noodles_bgzf as bgzf;
use noodles_csi::index::reference_sequence::bin::Chunk;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::Reader;
use crate::{Header, Record};

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
                State::Read(chunk_end) => match next_record(&mut ctx.reader, ctx.header).await? {
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

fn resolve_interval<B>(interval: B) -> (i32, i32)
where
    B: RangeBounds<i32>,
{
    match (interval.start_bound(), interval.end_bound()) {
        (Bound::Included(s), Bound::Included(e)) => (*s, *e),
        (Bound::Included(s), Bound::Unbounded) => (*s, i32::MAX),
        (Bound::Unbounded, Bound::Unbounded) => (1, i32::MAX),
        _ => todo!(),
    }
}

fn next_chunk(chunks: &[Chunk], i: &mut usize) -> Option<Chunk> {
    let chunk = chunks.get(*i).copied();
    *i += 1;
    chunk
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

fn intersects(
    record: &Record,
    reference_sequence_name: &str,
    interval_start: i32,
    interval_end: i32,
) -> io::Result<bool> {
    use crate::record::Chromosome;

    let name = match record.chromosome() {
        Chromosome::Name(s) => s.into(),
        Chromosome::Symbol(s) => s.to_string(),
    };

    let start = i32::from(record.position());
    let end = record
        .end()
        .map(i32::from)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(name == reference_sequence_name && in_interval(start, end, interval_start, interval_end))
}

fn in_interval(a_start: i32, a_end: i32, b_start: i32, b_end: i32) -> bool {
    a_start <= b_end && b_start <= a_end
}
