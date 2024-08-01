//! FASTQ writer.

mod builder;
mod record;

use std::io::{self, Write};

pub use self::builder::Builder;
use self::record::write_record;
use crate::Record;

const DEFAULT_DEFINITION_SEPARATOR: u8 = b' ';

/// A FASTQ writer.
pub struct Writer<W> {
    inner: W,
    definition_separator: u8,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a FASTQ writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    /// let writer = fastq::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self {
            inner,
            definition_separator: DEFAULT_DEFINITION_SEPARATOR,
        }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    /// let writer = fastq::io::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Writes a FASTQ record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fastq::{self as fastq, record::Definition};
    ///
    /// let mut writer = fastq::io::Writer::new(Vec::new());
    ///
    /// let record = fastq::Record::new(Definition::new("r0", ""), "ATCG", "NDLS");
    /// writer.write_record(&record)?;
    ///
    /// assert_eq!(writer.get_ref(), b"@r0\nATCG\n+\nNDLS\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        write_record(&mut self.inner, self.definition_separator, record)
    }
}
