use std::io::{self, Read, Seek};

use noodles_bgzf as bgzf;
use noodles_core::Region;

use super::{indexed_records::Record, Query};
use crate::BinningIndex;

/// An indexed reader.
pub struct IndexedReader<R, I> {
    inner: R,
    index: I,
}

impl<R, I> IndexedReader<bgzf::Reader<R>, I>
where
    R: Read,
    I: BinningIndex,
{
    /// Creates a indexed reader.
    pub fn new(inner: R, index: I) -> Self {
        Self {
            inner: bgzf::Reader::new(inner),
            index,
        }
    }

    /// Returns a reference to the underlying reader.
    pub fn get_ref(&self) -> &bgzf::Reader<R> {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    pub fn get_mut(&mut self) -> &mut bgzf::Reader<R> {
        &mut self.inner
    }

    /// Returns the underlying reader.
    pub fn into_inner(self) -> bgzf::Reader<R> {
        self.inner
    }

    /// Returns the associated index.
    pub fn index(&self) -> &I {
        &self.index
    }
}

impl<R, I> IndexedReader<bgzf::Reader<R>, I>
where
    R: Read + Seek,
    I: BinningIndex,
{
    /// Returns an iterator over records that intersects the given region.
    pub fn query<'r>(
        &'r mut self,
        region: &'r Region,
    ) -> io::Result<impl Iterator<Item = io::Result<Record>> + 'r> {
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

        Ok(Query::new(&mut self.inner, chunks)
            .indexed_records(header)
            .filter_by_region(region))
    }
}
