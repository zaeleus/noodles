use std::io::{self, Write};

use crate::fai::{Index, Record};

/// A FASTA index writer.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta::fai;
    /// let writer = fai::io::Writer::new(io::sink());
    /// let _inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta::fai;
    /// let mut writer = fai::io::Writer::new(io::sink());
    /// let _inner = writer.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta::fai;
    /// let writer = fai::io::Writer::new(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }
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
    /// let mut writer = fai::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Writes a FASTA index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use std::num::NonZero;
    ///
    /// use noodles_fasta::fai;
    ///
    /// let mut writer = fai::io::Writer::new(Vec::new());
    ///
    /// let line_base_count = const { NonZero::new(80).unwrap() };
    /// let line_width = const { NonZero::new(81).unwrap() };
    /// let index = fai::Index::from(vec![
    ///     fai::Record::new("sq0", 13, 5, line_base_count, line_width),
    /// ]);
    /// writer.write_index(&index)?;
    ///
    /// assert_eq!(writer.get_ref(), b"sq0\t13\t5\t80\t81\n");
    /// # Ok::<(), io::Error>(())
    pub fn write_index(&mut self, index: &Index) -> io::Result<()> {
        for record in index.as_ref() {
            write_record(&mut self.inner, record)?;
        }

        Ok(())
    }
}

pub(crate) fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(record.name())?;

    writeln!(
        writer,
        "\t{length}\t{offset}\t{line_base_count}\t{line_width}",
        length = record.length(),
        offset = record.position(),
        line_base_count = record.line_base_count(),
        line_width = record.line_width(),
    )
}
