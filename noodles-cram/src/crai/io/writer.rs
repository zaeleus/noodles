use std::io::{self, Write};

use flate2::write::GzEncoder;

use crate::crai::Record;

/// A CRAM index writer.
pub struct Writer<W>
where
    W: Write,
{
    inner: GzEncoder<W>,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a CRAM index writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram::crai;
    /// let writer = crai::io::Writer::new(io::sink());
    /// ```
    pub fn new(inner: W) -> Self {
        Self {
            inner: GzEncoder::new(inner, Default::default()),
        }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram::crai;
    /// let writer = crai::io::Writer::new(io::sink());
    /// let _inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &W {
        self.inner.get_ref()
    }

    /// Attempts to finish the output stream and returns the underlying writer.
    ///
    /// This is typically only manually called if the underlying stream is needed before the writer
    /// is dropped.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram::crai;
    /// let writer = crai::io::Writer::new(io::sink());
    /// let _inner = writer.finish()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn finish(self) -> io::Result<W> {
        self.inner.finish()
    }

    /// Writes a CRAM index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_core::Position;
    /// use noodles_cram::crai;
    ///
    /// let mut writer = crai::io::Writer::new(io::sink());
    ///
    /// let index = vec![crai::Record::new(
    ///     Some(0),
    ///     Position::new(10946),
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// )];
    ///
    /// writer.write_index(&index)?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn write_index(&mut self, index: &[Record]) -> io::Result<()> {
        write_index(&mut self.inner, index)
    }
}

fn write_index<W>(writer: &mut W, index: &[Record]) -> io::Result<()>
where
    W: Write,
{
    for record in index {
        write_record(writer, record)?;
    }

    Ok(())
}

fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: Write,
{
    writeln!(writer, "{record}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_core::Position;

        let index = vec![Record::new(
            Some(0),
            Position::new(10946),
            6765,
            17711,
            233,
            317811,
        )];

        let mut buf = Vec::new();
        write_index(&mut buf, &index)?;

        let expected = b"0\t10946\t6765\t17711\t233\t317811\n";

        assert_eq!(buf, expected);

        Ok(())
    }
}
