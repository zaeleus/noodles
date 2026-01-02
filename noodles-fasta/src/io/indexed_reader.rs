//! Indexed FASTA reader.

mod builder;
mod sequence;

pub use self::builder::Builder;

use std::io::{self, BufRead};

use noodles_core::Region;

use super::Reader;
use crate::{Record, fai};

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

    /// Returns the associated index.
    pub fn index(&self) -> &fai::Index {
        &self.index
    }
}

impl<R> IndexedReader<R>
where
    R: io::Read + io::Seek,
{
    /// Returns a record of the given region.
    ///
    /// This reads all bytes in one syscall and strips newlines in memory,
    /// which is significantly faster than line-by-line parsing for large sequences.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use noodles_core::Region;
    /// use noodles_fasta as fasta;
    ///
    /// let mut reader = fasta::io::indexed_reader::Builder::default()
    ///     .build_from_path("reference.fa")?;
    ///
    /// let region = "sq0:1-1000".parse()?;
    /// let record = reader.query(&region)?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query(&mut self, region: &Region) -> io::Result<Record> {
        use crate::record::{Definition, Sequence};
        use noodles_core::Position;

        let index_record = self
            .index
            .as_ref()
            .iter()
            .find(|r| r.name() == region.name())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("sequence not found: {}", region.name()),
                )
            })?;

        let interval = region.interval();
        let start = usize::from(interval.start().unwrap_or(Position::MIN));
        let end = interval
            .end()
            .map(usize::from)
            .unwrap_or(index_record.length() as usize);

        let start_base = (start - 1) as u64; // Convert 1-based to 0-based
        let len = (end - start + 1) as u64;

        let buf = sequence::read_sequence(self.inner.get_mut(), index_record, start_base, len)?;

        let definition = Definition::new(region.to_string(), None);
        let sequence = Sequence::from(buf);

        Ok(Record::new(definition, sequence))
    }
}
