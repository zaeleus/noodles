use std::io::{self, Write};

use byteorder::WriteBytesExt;
use noodles_bgzf as bgzf;

use super::MAGIC_NUMBER;

const MAJOR: u8 = 2;
const MINOR: u8 = 2;

/// A BCF writer.
pub struct Writer<W>
where
    W: Write,
{
    inner: bgzf::Writer<W>,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a BCF writer with a default compression level.
    ///
    /// The given stream is wrapped in a BGZF encoder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let writer = bcf::Writer::new(Vec::new());
    /// ```
    pub fn new(writer: W) -> Self {
        Self {
            inner: bgzf::Writer::new(writer),
        }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let writer = bcf::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        self.inner.get_ref()
    }

    /// Attempts to finish the output stream.
    ///
    /// This is typically only manually called if the underlying stream is needed before the writer
    /// is dropped.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    /// let mut writer = bcf::Writer::new(Vec::new());
    /// writer.try_finish()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn try_finish(&mut self) -> io::Result<()> {
        self.inner.try_finish()
    }

    /// Writes a BCF file format.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    /// let mut writer = bcf::Writer::new(Vec::new());
    /// writer.write_file_format()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_file_format(&mut self) -> io::Result<()> {
        self.inner.write_all(MAGIC_NUMBER)?;
        self.inner.write_u8(MAJOR)?;
        self.inner.write_u8(MINOR)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::Reader;

    use super::*;

    #[test]
    fn test_write_file_format() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());
        writer.write_file_format()?;
        writer.try_finish()?;

        let mut reader = Reader::new(writer.get_ref().as_slice());
        let file_format = reader.read_file_format()?;

        assert_eq!(file_format, (2, 2));

        Ok(())
    }
}
