use std::io::{self, BufRead, Seek};

use noodles_core::Region;

use crate::{fai, repository::Adapter, Reader, Record};

/// An indexed reader adapter.
pub struct IndexedReader<R> {
    reader: Reader<R>,
    index: fai::Index,
}

impl<R> IndexedReader<R>
where
    R: BufRead + Seek,
{
    /// Create a indexed reader adapter.
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
