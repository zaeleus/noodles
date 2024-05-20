//! FASTA writer.

mod builder;
mod record;

use std::io::{self, Write};

pub use self::builder::Builder;
use self::record::write_record;
use crate::Record;

/// A FASTA writer.
pub struct Writer<W> {
    inner: W,
    line_base_count: usize,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a FASTA writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// let writer = fasta::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Builder::default().build_with_writer(inner)
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// let writer = fasta::io::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Writes a FASTA record.
    ///
    /// By default, sequence lines are hard wrapped at 80 bases. This can be changed by using
    /// [`Builder::set_line_base_count`] when creating the writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta::{self as fasta, record::{Definition, Sequence}};
    ///
    /// let mut writer = fasta::io::Writer::new(Vec::new());
    ///
    /// let definition = Definition::new("sq0", None);
    /// let sequence = Sequence::from(b"ACGT".to_vec());
    /// let record = fasta::Record::new(definition, sequence);
    ///
    /// writer.write_record(&record)?;
    ///
    /// assert_eq!(writer.get_ref(), b">sq0\nACGT\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        write_record(&mut self.inner, record, self.line_base_count)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let writer = Writer::new(Vec::new());
        assert_eq!(writer.line_base_count, 80);
    }
}
