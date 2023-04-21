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
        write_record(&mut self.inner, record)
    }
}

fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(b"@")?;
    writer.write_all(record.name())?;

    if !record.description().is_empty() {
        writer.write_all(b" ")?;
        writer.write_all(record.description())?;
    }

    writer.write_all(b"\n")?;

    writer.write_all(record.sequence())?;
    writer.write_all(b"\n")?;

    writer.write_all(b"+\n")?;

    writer.write_all(record.quality_scores())?;
    writer.write_all(b"\n")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record() -> io::Result<()> {
        let mut record = Record::new("r0", "ACGT", "NDLS");

        let mut buf = Vec::new();
        write_record(&mut buf, &record)?;
        let expected = b"@r0\nACGT\n+\nNDLS\n";
        assert_eq!(buf, expected);

        record.description_mut().extend_from_slice(b"LN:4");

        buf.clear();
        write_record(&mut buf, &record)?;
        let expected = b"@r0 LN:4\nACGT\n+\nNDLS\n";
        assert_eq!(buf, expected);

        Ok(())
    }
}
