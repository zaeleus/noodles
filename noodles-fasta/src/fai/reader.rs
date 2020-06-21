use std::io::{self, BufRead};

/// A FASTA index reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    /// Creates a FASTA index reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// let data = b"sq0\t13\t5\t80\t81\nsq1\t21\t19\t80\t81\n";
    /// let mut reader = fai::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a raw FASTA index record.
    ///
    /// The given buffer will not include the trailing newline. It can subsequently be parsed as a
    /// [`fai::Record`].
    ///
    /// The position of the stream is expected to be at the start or at the start of another
    /// record.
    ///
    /// If successful, this returns the number of bytes read from the stream. If the number of
    /// bytes read is 0, the stream reached EOF.
    ///
    /// [`fai::Record`]: struct.Record.html
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta::fai;
    ///
    /// let data = b"sq0\t13\t5\t80\t81\nsq1\t21\t19\t80\t81\n";
    /// let mut reader = fai::Reader::new(&data[..]);
    ///
    /// let mut buf = String::new();
    /// reader.read_record(&mut buf)?;
    ///
    /// assert_eq!(buf, "sq0\t13\t5\t80\t81");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_record(&mut self, buf: &mut String) -> io::Result<usize> {
        let result = self.inner.read_line(buf);
        buf.pop();
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_record() -> io::Result<()> {
        let data = b"\
sq0\t10946\t4\t80\t81
sq1\t17711\t10954\t80\t81
";

        let mut reader = Reader::new(&data[..]);

        let mut buf = String::new();
        reader.read_record(&mut buf)?;
        assert_eq!(buf, "sq0\t10946\t4\t80\t81");

        buf.clear();
        reader.read_record(&mut buf)?;
        assert_eq!(buf, "sq1\t17711\t10954\t80\t81");

        let bytes_read = reader.read_record(&mut buf)?;
        assert_eq!(bytes_read, 0);

        Ok(())
    }
}
