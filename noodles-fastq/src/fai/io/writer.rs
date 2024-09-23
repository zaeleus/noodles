use std::io::{self, Write};

use crate::fai::Record;

/// A FASTQ index writer.
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
    /// use noodles_fastq::fai;
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
    /// use noodles_fastq::fai;
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
    /// use noodles_fastq::fai;
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
    /// Creates a FASTQ index writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq::fai;
    /// let mut writer = fai::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Writes a FASTQ index record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fastq::fai;
    ///
    /// let mut writer = fai::io::Writer::new(Vec::new());
    ///
    /// let record = fai::Record::new("r0", 4, 4, 4, 5, 11);
    /// writer.write_record(&record)?;
    ///
    /// assert_eq!(writer.get_ref(), b"r0\t4\t4\t4\t5\t11\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        writeln!(
            self.inner,
            "{name}\t{length}\t{sequence_offset}\t{line_bases}\t{line_width}\t{quality_scores_offset}",
            name = record.name(),
            length = record.length(),
            sequence_offset = record.sequence_offset(),
            line_bases = record.line_bases(),
            line_width = record.line_width(),
            quality_scores_offset = record.quality_scores_offset(),
        )
    }
}
