use std::vec;

use futures::{Stream, stream};
use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::Reader;
use crate::Record;

enum State {
    Seek,
    Read(bgzf::VirtualPosition),
    Done,
}

/// An async reader over records of an async BCF reader that intersects a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'r, R>
where
    R: AsyncRead + AsyncSeek,
{
    inner: &'r mut Reader<bgzf::r#async::io::Reader<R>>,
    chunks: vec::IntoIter<Chunk>,

    reference_sequence_id: usize,
    interval: Interval,

    state: State,
}

impl<'r, R> Query<'r, R>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    pub(super) fn new(
        inner: &'r mut Reader<bgzf::r#async::io::Reader<R>>,
        chunks: Vec<Chunk>,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            inner,
            chunks: chunks.into_iter(),
            reference_sequence_id,
            interval,
            state: State::Seek,
        }
    }

    pub async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
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

                        if intersects(record, self.reference_sequence_id, self.interval)? {
                            return Ok(n);
                        }
                    }
                },
                State::Done => return Ok(0),
            }
        }
    }

    pub fn records(self) -> impl Stream<Item = io::Result<Record>> {
        Box::pin(stream::try_unfold(self, |mut reader| async {
            let mut record = Record::default();

            match reader.read_record(&mut record).await? {
                0 => Ok(None),
                _ => Ok(Some((record, reader))),
            }
        }))
    }
}

fn intersects(
    record: &Record,
    chromosome_id: usize,
    region_interval: Interval,
) -> io::Result<bool> {
    let id = record.reference_sequence_id()?;

    let Some(start) = record.variant_start().transpose()? else {
        return Ok(false);
    };

    let end = record.end()?;

    let record_interval = Interval::from(start..=end);

    Ok(id == chromosome_id && record_interval.intersects(region_interval))
}
