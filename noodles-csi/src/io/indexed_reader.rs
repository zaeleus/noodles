use std::io::{self, Read, Seek};

use noodles_bgzf as bgzf;
use noodles_core::Region;

use super::{indexed_records::Record, Query};
use crate::{BinningIndex, Index};

/// An indexed reader.
pub struct IndexedReader<R> {
    inner: R,
    index: Index,
}

impl<R> IndexedReader<bgzf::Reader<R>>
where
    R: Read,
{
    /// Creates a indexed reader.
    pub fn new(inner: R, index: Index) -> Self {
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
    pub fn index(&self) -> &Index {
        &self.index
    }
}

impl<R> IndexedReader<bgzf::Reader<R>>
where
    R: Read + Seek,
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
