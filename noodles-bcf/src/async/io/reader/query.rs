use futures::{Stream, stream};
use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi as csi;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::Reader;
use crate::Record;

/// An async reader over records of an async BCF reader that intersects a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'r, R>
where
    R: AsyncRead + AsyncSeek,
{
    reader: Reader<csi::r#async::io::Query<'r, R>>,
    reference_sequence_id: usize,
    interval: Interval,
}

impl<'r, R> Query<'r, R>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    pub(super) fn new(
        inner: &'r mut bgzf::r#async::io::Reader<R>,
        chunks: Vec<Chunk>,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            reader: Reader::from(csi::r#async::io::Query::new(inner, chunks)),
            reference_sequence_id,
            interval,
        }
    }

    pub async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        loop {
            match self.reader.read_record(record).await? {
                0 => return Ok(0),
                n => {
                    if intersects(record, self.reference_sequence_id, self.interval)? {
                        return Ok(n);
                    }
                }
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
