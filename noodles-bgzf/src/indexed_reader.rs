use std::io::{self, Read, Seek, SeekFrom};

use super::{gzi, Reader};

/// An indexed BGZF reader.
pub struct IndexedReader<R> {
    inner: Reader<R>,
    index: gzi::Index,
}

impl<R> IndexedReader<R>
where
    R: Read,
{
    /// Creates an indexed BGZF reader.
    pub fn new(inner: R, index: gzi::Index) -> Self {
        Self {
            inner: Reader::new(inner),
            index,
        }
    }
}

impl<R> Read for IndexedReader<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buf)
    }
}

impl<R> Seek for IndexedReader<R>
where
    R: Read + Seek,
{
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        match pos {
            SeekFrom::Start(pos) => self.inner.seek_by_uncompressed_position(&self.index, pos),
            _ => unimplemented!(),
        }
    }
}
