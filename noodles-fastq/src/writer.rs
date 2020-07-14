use std::io::{self, Write};

use super::Record;

/// A FASTQ writer.
pub struct Writer<W> {
    inner: W,
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
    /// let writer = fastq::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastq as fastq;
    /// let writer = fastq::Writer::new(Vec::new());
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
    /// use noodles_fastq as fastq;
    ///
    /// let mut writer = fastq::Writer::new(Vec::new());
    ///
    /// let record = fastq::Record::new("r0", "ATCG", "NDLS");
    /// writer.write_record(&record)?;
    ///
    /// assert_eq!(writer.get_ref(), b"@r0\nATCG\n+\nNDLS\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        self.inner.write_all(b"@")?;
        self.inner.write_all(record.read_name())?;
        self.inner.write_all(b"\n")?;
        self.inner.write_all(record.sequence())?;
        self.inner.write_all(b"\n+\n")?;
        self.inner.write_all(record.quality_scores())?;
        self.inner.write_all(b"\n")?;

        Ok(())
    }
}
