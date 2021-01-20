use std::io::{self, Write};

use super::Record;

/// A FASTA index writer.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a FASTA index writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// let mut writer = fai::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// let mut writer = fai::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Writes a FASTA index record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta::fai;
    ///
    /// let mut writer = fai::Writer::new(Vec::new());
    ///
    /// let record = fai::Record::new(String::from("sq0"), 13, 5, 80, 81);
    /// writer.write_record(&record)?;
    ///
    /// assert_eq!(writer.get_ref(), b"sq0\t13\t5\t80\t81\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        writeln!(
            self.inner,
            "{reference_sequence_name}\t{len}\t{offset}\t{line_bases}\t{line_width}",
            reference_sequence_name = record.reference_sequence_name(),
            len = record.len(),
            offset = record.offset(),
            line_bases = record.line_bases(),
            line_width = record.line_width(),
        )
    }
}
