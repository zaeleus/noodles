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

    /// Writes a FASTA index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta::fai;
    ///
    /// let mut writer = fai::Writer::new(Vec::new());
    ///
    /// let index = vec![fai::Record::new(String::from("sq0"), 13, 5, 80, 81)];
    /// writer.write_index(&index)?;
    ///
    /// assert_eq!(writer.get_ref(), b"sq0\t13\t5\t80\t81\n");
    /// # Ok::<(), io::Error>(())
    pub fn write_index(&mut self, index: &[Record]) -> io::Result<()> {
        for record in index {
            write_record(&mut self.inner, record)?;
        }

        Ok(())
    }
}

fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: Write,
{
    writeln!(
        writer,
        "{name}\t{len}\t{offset}\t{line_bases}\t{line_width}",
        name = record.name(),
        len = record.len(),
        offset = record.offset(),
        line_bases = record.line_bases(),
        line_width = record.line_width(),
    )
}
