use futures::{Stream, stream};
use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi as csi;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::Reader;
use crate::{Record, io::reader::query::intersects};

/// An async reader over records of an async BAM reader that intersects a given region.
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
        reader: &'r mut bgzf::r#async::io::Reader<R>,
        chunks: Vec<Chunk>,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            reader: Reader::from(csi::r#async::io::Query::new(reader, chunks)),
            reference_sequence_id,
            interval,
        }
    }

    /// Reads a record.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use noodles_bam::{self as bam, bai};
    /// use noodles_core::Region;
    /// use noodles_sam as sam;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam").await.map(bam::r#async::io::Reader::new)?;
    /// let header = reader.read_header().await?;
    ///
    /// let index = bai::r#async::fs::read("sample.bam.bai").await?;
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(&header, &index, &region)?;
    ///
    /// let mut record = bam::Record::default();
    ///
    /// while query.read_record(&mut record).await? != 0 {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
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
