//! Indexed BGZF reader.

mod builder;

pub use self::builder::Builder;

use std::io::{self, BufRead, Read, Seek, SeekFrom};

use super::{gzi, Reader, VirtualPosition};

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

    /// Returns a reference to the underlying reader.
    pub fn get_ref(&self) -> &R {
        self.inner.get_ref()
    }

    /// Returns a mutable reference to the underlying reader.
    pub fn get_mut(&mut self) -> &mut R {
        self.inner.get_mut()
    }

    /// Unwraps and returns the underlying writer.
    pub fn into_inner(self) -> R {
        self.inner.into_inner()
    }

    /// Returns the current position of the stream.
    pub fn position(&self) -> u64 {
        self.inner.position()
    }

    /// Returns the current virtual position of the stream.
    pub fn virtual_position(&self) -> VirtualPosition {
        self.inner.virtual_position()
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

impl<R> BufRead for IndexedReader<R>
where
    R: Read,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.inner.fill_buf()
    }

    fn consume(&mut self, amt: usize) {
        self.inner.consume(amt)
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
