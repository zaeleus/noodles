//! FASTA I/O.

pub mod indexed_reader;
mod indexer;
pub mod reader;
pub mod writer;

use std::io::{self, BufRead, Read, Seek, SeekFrom};

use noodles_bgzf as bgzf;

pub use self::{indexed_reader::IndexedReader, indexer::Indexer, reader::Reader, writer::Writer};

#[deprecated(since = "0.48.0", note = "Use `fasta::fs::index` instead.")]
pub use super::fs::index;

/// A buffered FASTA reader.
pub enum BufReader<R> {
    /// bgzip-compressed.
    Bgzf(bgzf::IndexedReader<R>),
    /// Uncompressed.
    Uncompressed(std::io::BufReader<R>),
}

impl<R> Read for BufReader<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Self::Bgzf(reader) => reader.read(buf),
            Self::Uncompressed(reader) => reader.read(buf),
        }
    }
}

impl<R> BufRead for BufReader<R>
where
    R: Read,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            Self::Bgzf(reader) => reader.fill_buf(),
            Self::Uncompressed(reader) => reader.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            Self::Bgzf(reader) => reader.consume(amt),
            Self::Uncompressed(reader) => reader.consume(amt),
        }
    }
}

impl<R> Seek for BufReader<R>
where
    R: Read + Seek,
{
    fn seek(&mut self, pos: SeekFrom) -> std::io::Result<u64> {
        match self {
            Self::Bgzf(reader) => reader.seek(pos),
            Self::Uncompressed(reader) => reader.seek(pos),
        }
    }
}
