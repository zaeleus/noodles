use std::io::{self, Read, Seek};

use noodles_bgzf as bgzf;
use noodles_core::Region;

use super::{Query, indexed_records::Record};
use crate::BinningIndex;

/// An indexed reader.
pub struct IndexedReader<R, I> {
    inner: R,
    index: I,
}

impl<R, I> IndexedReader<bgzf::io::Reader<R>, I>
where
    R: Read,
    I: BinningIndex,
{
    /// Creates an indexed reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    /// let index = csi::Index::default();
    /// let reader = csi::io::IndexedReader::new(io::empty(), index);
    /// ```
    pub fn new(inner: R, index: I) -> Self {
        Self {
            inner: bgzf::io::Reader::new(inner),
            index,
        }
    }

    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    ///
    /// let index = csi::Index::default();
    /// let reader = csi::io::IndexedReader::new(io::empty(), index);
    /// let inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &bgzf::io::Reader<R> {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    ///
    /// let index = csi::Index::default();
    /// let mut reader = csi::io::IndexedReader::new(io::empty(), index);
    /// let inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut bgzf::io::Reader<R> {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    ///
    /// let index = csi::Index::default();
    /// let reader = csi::io::IndexedReader::new(io::empty(), index);
    /// let inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> bgzf::io::Reader<R> {
        self.inner
    }

    /// Returns the associated index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_csi as csi;
    ///
    /// let index = csi::Index::default();
    /// let reader = csi::io::IndexedReader::new(io::empty(), index);
    /// let index = reader.index();
    /// ```
    pub fn index(&self) -> &I {
        &self.index
    }
}

impl<R, I> IndexedReader<bgzf::io::Reader<R>, I>
where
    R: Read + Seek,
    I: BinningIndex,
{
    /// Returns an iterator over records that intersects the given region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_csi as csi;
    ///
    /// let index = csi::fs::read("sample.bed.gz.csi")?;
    /// let mut reader = File::open("sample.bed.gz")
    ///     .map(|f| csi::io::IndexedReader::new(f, index))?;
    ///
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&region)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
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
