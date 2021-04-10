use std::io::{self, Read};

use noodles_bgzf as bgzf;

static MAGIC_NUMBER: &[u8] = b"BCF";

/// A BCF reader.
///
/// The BCF format is comprised of two parts: 1) a VCF header and 2) a list of records.
pub struct Reader<R> {
    inner: bgzf::Reader<R>,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates a BCF reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bcf as bcf;
    /// let mut reader = File::open("sample.bcf").map(bcf::Reader::new)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn new(reader: R) -> Self {
        Self {
            inner: bgzf::Reader::new(reader),
        }
    }

    /// Reads the BCF file format.
    ///
    /// The BCF magic number is also checked.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// This returns the major and minor format versions as a tuple.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bcf as bcf;
    /// let mut reader = File::open("sample.bcf").map(bcf::Reader::new)?;
    /// let (major, minor) = reader.read_file_format()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_file_format(&mut self) -> io::Result<(u8, u8)> {
        let mut buf = [0; 5];

        self.inner.read_exact(&mut buf)?;

        let magic = &buf[..3];

        if magic != MAGIC_NUMBER {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid file format",
            ));
        }

        let major = buf[3];
        let minor = buf[4];

        Ok((major, minor))
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    fn compress(data: &[u8]) -> io::Result<Vec<u8>> {
        let mut writer = bgzf::Writer::new(Vec::new());
        writer.write_all(data)?;
        writer.finish()
    }

    #[test]
    fn test_read_file_format() -> io::Result<()> {
        let data = compress(b"BCF\x02\x01")?;
        let mut reader = Reader::new(&data[..]);

        let (major, minor) = reader.read_file_format()?;

        assert_eq!(major, 2);
        assert_eq!(minor, 1);

        Ok(())
    }

    #[test]
    fn test_read_file_format_with_an_invalid_magic_number() -> io::Result<()> {
        let data = compress(b"BAM\x02\x01")?;
        let mut reader = Reader::new(&data[..]);
        assert!(reader.read_file_format().is_err());
        Ok(())
    }
}
