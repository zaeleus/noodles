mod builder;

pub use self::builder::Builder;

use std::{
    fs::File,
    io::{self, BufRead, BufReader, Seek},
};

use noodles_core::Region;

use crate::{fai, repository::Adapter, Reader, Record};

/// An indexed reader adapter.
pub struct IndexedReader<R> {
    reader: Reader<R>,
    index: fai::Index,
}

impl IndexedReader<BufReader<File>> {
    /// Creates an indexed reader adapter builder for paths on a filesystem.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::repository::adapters::IndexedReader;
    /// let builder = IndexedReader::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }
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
    /// let reader = fasta::Reader::new(io::empty());
    /// let index = fai::Index::default();
    /// let adapter = IndexedReader::new(reader, index);
    /// ```
    pub fn new(reader: Reader<R>, index: fai::Index) -> Self {
        Self { reader, index }
    }
}

impl<R> Adapter for IndexedReader<R>
where
    R: BufRead + Seek,
{
    fn get(&mut self, name: &str) -> Option<io::Result<Record>> {
        let region = Region::new(name, ..);
        Some(self.reader.query(&self.index, &region))
    }
}
