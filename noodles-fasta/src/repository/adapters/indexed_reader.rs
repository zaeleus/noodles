use std::io::{self, BufRead, Seek};

use noodles_core::Region;

use crate::{repository::Adapter, Record};

/// An indexed reader adapter.
pub struct IndexedReader<R> {
    reader: crate::io::IndexedReader<R>,
}

impl<R> IndexedReader<R>
where
    R: BufRead + Seek,
{
    /// Creates an indexed reader adapter.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta::{self as fasta, fai, repository::adapters::IndexedReader};
    /// let index = fai::Index::default();
    /// let reader = fasta::io::IndexedReader::new(io::empty(), index);
    /// let adapter = IndexedReader::new(reader);
    /// ```
    pub fn new(reader: crate::io::IndexedReader<R>) -> Self {
        Self { reader }
    }
}

impl<R> Adapter for IndexedReader<R>
where
    R: BufRead + Seek + Send + Sync,
{
    fn get(&mut self, name: &[u8]) -> Option<io::Result<Record>> {
        let region = Region::new(name, ..);
        Some(self.reader.query(&region))
    }
}
