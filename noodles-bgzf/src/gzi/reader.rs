use std::io;
use std::io::Read;
use super::Index;

/// A GZI index reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
    where
        R: Read,
{
    /// Creates a GZI index reader.
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a GZI index.
    ///
    /// The position of the [`Read`](std::io::Read) stream is expected to be at the start.
    pub fn read_index(&mut self) -> io::Result<Index> {
        todo!()
    }
}