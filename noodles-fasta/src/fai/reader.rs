use std::io::{self, BufRead};

use super::Record;

/// A FASTA index reader.
pub struct Reader<R> {
    inner: R,
    line_buf: String,
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
        Self {
            inner,
            line_buf: String::new(),
        }
    }

    /// Reads a FASTA index record.
    ///
    /// The position of the stream is expected to be at the start or at the start of another
    /// record.
    ///
    /// If successful, this returns the number of bytes read from the stream. If the number of
    /// bytes read is 0, the stream reached EOF.
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
    /// let mut record = fai::Record::default();
    /// reader.read_record(&mut record)?;
    ///
    /// assert_eq!(record.name(), "sq0");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        self.line_buf.clear();

        match self.inner.read_line(&mut self.line_buf) {
            Ok(0) => Ok(0),
            Ok(n) => {
                self.line_buf.pop();

                *record = self
                    .line_buf
                    .parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                Ok(n)
            }
            Err(e) => Err(e),
        }
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
        let mut record = Record::default();

        let bytes_read = reader.read_record(&mut record)?;
        assert_eq!(bytes_read, 18);

        let bytes_read = reader.read_record(&mut record)?;
        assert_eq!(bytes_read, 22);

        let bytes_read = reader.read_record(&mut record)?;
        assert_eq!(bytes_read, 0);

        Ok(())
    }
}
