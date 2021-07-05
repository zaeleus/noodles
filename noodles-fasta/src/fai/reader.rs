use std::io::{self, BufRead};

use super::Index;

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

    /// Reads a FASTA index.
    ///
    /// The position of the stream is expected to be at the start or at the start of a record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta::fai;
    ///
    /// let data = b"sq0\t13\t5\t80\t81\nsq1\t21\t19\t80\t81\n";
    /// let mut reader = fai::Reader::new(&data[..]);
    /// let index = reader.read_index()?;
    ///
    /// assert_eq!(index, vec![
    ///     fai::Record::new(String::from("sq0"), 13, 5, 80, 81),
    ///     fai::Record::new(String::from("sq1"), 21, 19, 80, 81),
    /// ]);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_index(&mut self) -> io::Result<Index> {
        let mut buf = String::new();
        let mut index = Vec::new();

        loop {
            buf.clear();

            match read_line(&mut self.inner, &mut buf) {
                Ok(0) => break,
                Ok(_) => {
                    let record = buf
                        .parse()
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                    index.push(record);
                }
                Err(e) => return Err(e),
            }
        }

        Ok(index)
    }
}

pub fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    let result = reader.read_line(buf);
    buf.pop();
    result
}
