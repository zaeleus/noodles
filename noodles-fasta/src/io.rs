//! FASTA I/O.

pub mod indexed_reader;
mod indexer;
pub mod reader;
pub mod writer;

use std::{
    fs::File,
    io::{self, BufRead, Read, Seek, SeekFrom},
    path::Path,
};

use noodles_bgzf as bgzf;

pub use self::{indexed_reader::IndexedReader, indexer::Indexer, reader::Reader, writer::Writer};
use super::fai;

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

/// Indexes a FASTA file.
///
/// # Examples
///
/// ```no_run
/// # use std::io;
/// use noodles_fasta as fasta;
/// let index = fasta::io::index("reference.fa")?;
/// # Ok::<(), io::Error>(())
/// ```
pub fn index<P>(src: P) -> io::Result<fai::Index>
where
    P: AsRef<Path>,
{
    let mut indexer = File::open(src).map(io::BufReader::new).map(Indexer::new)?;
    let mut index = Vec::new();

    while let Some(i) = indexer.index_record()? {
        index.push(i);
    }

    Ok(index)
}
