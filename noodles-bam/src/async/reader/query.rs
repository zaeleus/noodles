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

struct Context {
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

    let context = Context {
        chunks,
        i: 0,

        reference_sequence_id,
        start,
        end,

        state: State::Seek,
    };

    Box::pin(stream::unfold(
        (reader, context),
        |(mut reader, mut context)| async {
            loop {
                match context.state {
                    State::Seek => {
                        context.state = match next_chunk(&mut context) {
                            Some(chunk) => {
                                if let Err(e) = reader.seek(chunk.start()).await {
                                    return Some((Err(e), (reader, context)));
                                }

                                State::Read(chunk.end())
                            }
                            None => State::Done,
                        };
                    }
                    State::Read(chunk_end) => match next_record(&mut reader).await {
                        Some(Ok(record)) => {
                            if reader.virtual_position() >= chunk_end {
                                context.state = State::Seek;
                            }

                            if let Some(reference_sequence_id) = record.reference_sequence_id() {
                                let reference_sequence_id =
                                    i32::from(reference_sequence_id) as usize;

                                let start =
                                    record.position().map(i32::from).expect("missing position");
                                let len = match record.cigar().reference_len() {
                                    Ok(len) => len as i32,
                                    Err(e) => return Some((Err(e), (reader, context))),
                                };
                                let end = start + len - 1;

                                if reference_sequence_id == context.reference_sequence_id
                                    && in_interval(start, end, context.start, context.end)
                                {
                                    return Some((Ok(record), (reader, context)));
                                }
                            }
                        }
                        Some(Err(e)) => {
                            return Some((Err(e), (reader, context)));
                        }
                        None => {
                            context.state = State::Seek;
                        }
                    },
                    State::Done => {
                        return None;
                    }
                }
            }
        },
    ))
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

fn next_chunk(context: &mut Context) -> Option<Chunk> {
    let chunk = context.chunks.get(context.i).copied();
    context.i += 1;
    chunk
}

async fn next_record<R>(reader: &mut Reader<R>) -> Option<io::Result<Record>>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    let mut record = Record::default();

    match reader.read_record(&mut record).await {
        Ok(0) => None,
        Ok(_) => Some(Ok(record)),
        Err(e) => Some(Err(e)),
    }
}

fn in_interval(a_start: i32, a_end: i32, b_start: i32, b_end: i32) -> bool {
    a_start <= b_end && b_start <= a_end
}
