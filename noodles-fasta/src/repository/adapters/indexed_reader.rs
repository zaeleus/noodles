use std::io::{self, BufRead, Seek};

use noodles_core::Region;

use crate::{repository::Adapter, Reader, Record};

/// An indexed reader adapter.
pub struct IndexedReader<R> {
    reader: Reader<R>,
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
    /// use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
    /// let reader = fasta::Reader::new(io::empty());
    /// let adapter = IndexedReader::new(reader);
    /// ```
    pub fn new(reader: Reader<R>) -> Self {
        Self { reader }
    }
}

impl<R> Adapter for IndexedReader<R>
where
    R: BufRead + Seek,
{
    fn get(&mut self, name: &str) -> Option<io::Result<Record>> {
        let region = Region::new(name, ..);
        Some(self.reader.query(&region))
    }
}
