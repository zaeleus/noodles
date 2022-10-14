//! Indexed FASTA reader.

mod builder;

pub use self::builder::Builder;

use std::io::{self, BufRead, Seek};

use noodles_core::Region;

use super::{fai, Reader, Record};

/// An indexed FASTA reader.
pub struct IndexedReader<R> {
    inner: Reader<R>,
    index: fai::Index,
}

impl<R> IndexedReader<R>
where
    R: BufRead,
{
    /// Creates a new indexed FASTA reader.
    pub fn new(inner: R, index: fai::Index) -> Self {
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

    /// Returns the underlying reader.
    pub fn into_inner(self) -> R {
        self.inner.into_inner()
    }

    /// Reads a raw definition line.
    pub fn read_definition(&mut self, buf: &mut String) -> io::Result<usize> {
        self.inner.read_definition(buf)
    }

    /// Reads a sequence.
    pub fn read_sequence(&mut self, buf: &mut Vec<u8>) -> io::Result<usize> {
        self.inner.read_sequence(buf)
    }
}

impl<R> IndexedReader<R>
where
    R: BufRead + Seek,
{
    /// Returns a record of the given region.
    pub fn query(&mut self, region: &Region) -> io::Result<Record> {
        self.inner.query(&self.index, region)
    }
}
