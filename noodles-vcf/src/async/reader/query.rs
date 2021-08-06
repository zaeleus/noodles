use std::ops::{Bound, RangeBounds};

use futures::{stream, Stream};
use noodles_bgzf as bgzf;
use noodles_csi::index::reference_sequence::bin::Chunk;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::Reader;
use crate::Record;

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    Done,
}

struct Context<'a, R>
where
    R: AsyncRead + AsyncSeek,
{
    reader: &'a mut Reader<bgzf::AsyncReader<R>>,

    chunks: Vec<Chunk>,
    i: usize,

    reference_sequence_name: String,
    start: i32,
    end: i32,

    state: State,
}

pub fn query<R, B>(
    reader: &mut Reader<bgzf::AsyncReader<R>>,
    chunks: Vec<Chunk>,
    reference_sequence_name: String,
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

        reference_sequence_name,
        start,
        end,

        state: State::Seek,
    };

    Box::pin(stream::unfold(ctx, |mut ctx| async {
        loop {
            match ctx.state {
                State::Seek => {
                    ctx.state = match next_chunk(&ctx.chunks, &mut ctx.i) {
                        Some(chunk) => {
                            if let Err(e) = ctx.reader.seek(chunk.start()).await {
                                return Some((Err(e), ctx));
                            }

                            State::Read(chunk.end())
                        }
                        None => State::Done,
                    };
                }
                State::Read(chunk_end) => match next_record(&mut ctx.reader).await {
                    Some(Ok(record)) => {
                        if ctx.reader.virtual_position() >= chunk_end {
                            ctx.state = State::Seek;
                        }

                        match intersects(&record, &ctx.reference_sequence_name, ctx.start, ctx.end)
                        {
                            Ok(true) => return Some((Ok(record), ctx)),
                            Ok(false) => {}
                            Err(e) => return Some((Err(e), ctx)),
                        }
                    }
                    Some(Err(e)) => return Some((Err(e), ctx)),
                    None => ctx.state = State::Seek,
                },
                State::Done => return None,
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

async fn next_record<R>(reader: &mut Reader<bgzf::AsyncReader<R>>) -> Option<io::Result<Record>>
where
    R: AsyncRead + Unpin,
{
    let mut buf = String::new();

    match reader.read_record(&mut buf).await {
        Ok(0) => None,
        Ok(_) => {
            let result = buf
                .parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e));

            Some(result)
        }
        Err(e) => Some(Err(e)),
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
