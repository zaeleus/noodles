use futures::Stream;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use tokio::io::{self, AsyncRead, AsyncSeek};

use super::{Query, filtered_indexed_records};
use crate::{BinningIndex, io::indexed_records::Record};

/// An async indexed reader.
pub struct IndexedReader<R, I> {
    inner: R,
    index: I,
}

impl<R, I> IndexedReader<bgzf::r#async::io::Reader<R>, I>
where
    R: AsyncRead,
    I: BinningIndex,
{
    /// Creates an async indexed reader.
    pub fn new(inner: R, index: I) -> Self {
        Self {
            inner: bgzf::r#async::io::Reader::new(inner),
            index,
        }
    }

    /// Returns the associated index.
    pub fn index(&self) -> &I {
        &self.index
    }
}

impl<R, I> IndexedReader<bgzf::r#async::io::Reader<R>, I>
where
    R: AsyncRead + AsyncSeek + Unpin,
    I: BinningIndex,
{
    /// Returns a stream over records that intersects the given region.
    pub fn query(&mut self, region: &Region) -> io::Result<impl Stream<Item = io::Result<Record>>> {
        let header = self
            .index
            .header()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing index header"))?;

        let reference_sequence_id = header
            .reference_sequence_names()
            .get_index_of(region.name())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "missing reference sequence name",
                )
            })?;

        let chunks = self.index.query(reference_sequence_id, region.interval())?;

        Ok(filtered_indexed_records(
            Query::new(&mut self.inner, chunks),
            header,
            region,
        ))
    }
}
