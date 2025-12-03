use futures::{Stream, stream};
use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi as csi;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::{Header, Reader};
use crate::{Record, io::reader::query::intersects};

/// An async reader over records of an async SAM reader that intersects a given a region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'r, 'h: 'r, R>
where
    R: AsyncRead + AsyncSeek,
{
    reader: Reader<csi::r#async::io::Query<'r, R>>,
    header: &'h Header,
    reference_sequence_id: usize,
    interval: Interval,
}

impl<'r, 'h: 'r, R> Query<'r, 'h, R>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    pub(super) fn new(
        reader: &'r mut bgzf::r#async::io::Reader<R>,
        chunks: Vec<Chunk>,
        header: &'h Header,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            reader: Reader::new(csi::r#async::io::Query::new(reader, chunks)),
            header,
            reference_sequence_id,
            interval,
        }
    }

    async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        loop {
            match self.reader.read_record(record).await? {
                0 => return Ok(0),
                n => {
                    if intersects(
                        self.header,
                        record,
                        self.reference_sequence_id,
                        self.interval,
                    )? {
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
