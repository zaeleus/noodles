use tokio::io::{self, AsyncBufRead};

use crate::fai::Index;

/// An async FASTA index reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: AsyncBufRead + Unpin,
{
    /// Creates an async FASTA index reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// let data = [];
    /// let mut reader = fai::AsyncReader::new(&data[..]);
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
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_fasta::fai;
    ///
    /// let data = b"sq0\t13\t5\t80\t81\nsq1\t21\t19\t80\t81\n";
    /// let mut reader = fai::AsyncReader::new(&data[..]);
    ///
    /// let index = reader.read_index().await?;
    ///
    /// assert_eq!(index, vec![
    ///     fai::Record::new(String::from("sq0"), 13, 5, 80, 81),
    ///     fai::Record::new(String::from("sq1"), 21, 19, 80, 81),
    /// ]);
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_index(&mut self) -> io::Result<Index> {
        use crate::r#async::reader::read_line;

        let mut buf = String::new();
        let mut index = Vec::new();

        loop {
            buf.clear();

            match read_line(&mut self.inner, &mut buf).await? {
                0 => break,
                _ => {
                    let record = buf
                        .parse()
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                    index.push(record);
                }
            }
        }

        Ok(index)
    }
}
