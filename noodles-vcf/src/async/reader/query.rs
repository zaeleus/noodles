use std::vec;

use futures::{stream, Stream};
use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::index::reference_sequence::bin::Chunk;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::Reader;
use crate::{reader::query::intersects, Header, Record};

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

    chunks: vec::IntoIter<Chunk>,

    reference_sequence_name: String,
    interval: Interval,

    state: State,

    header: &'r Header,
}

pub fn query<'r, R>(
    reader: &'r mut Reader<bgzf::AsyncReader<R>>,
    chunks: Vec<Chunk>,
    reference_sequence_name: String,
    interval: Interval,
    header: &'r Header,
) -> impl Stream<Item = io::Result<Record>> + 'r
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    let ctx = Context {
        reader,

        chunks: chunks.into_iter(),

        reference_sequence_name,
        interval,

        state: State::Seek,

        header,
    };

    Box::pin(stream::try_unfold(ctx, |mut ctx| async {
        loop {
            match ctx.state {
                State::Seek => {
                    ctx.state = match ctx.chunks.next() {
                        Some(chunk) => {
                            ctx.reader.seek_virtual_position(chunk.start()).await?;
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

                        if intersects(&record, &ctx.reference_sequence_name, ctx.interval)? {
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
