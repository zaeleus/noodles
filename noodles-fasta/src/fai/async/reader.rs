use tokio::io::{self, AsyncBufRead, AsyncBufReadExt};

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

async fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    const LINE_FEED: char = '\n';
    const CARRIAGE_RETURN: char = '\r';

    match reader.read_line(buf).await? {
        0 => Ok(0),
        n => {
            if buf.ends_with(LINE_FEED) {
                buf.pop();

                if buf.ends_with(CARRIAGE_RETURN) {
                    buf.pop();
                }
            }

            Ok(n)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_line() -> io::Result<()> {
        async fn t(buf: &mut String, mut data: &[u8], expected: &str) -> io::Result<()> {
            buf.clear();
            read_line(&mut data, buf).await?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = String::new();

        t(&mut buf, b"noodles\n", "noodles").await?;
        t(&mut buf, b"noodles\r\n", "noodles").await?;
        t(&mut buf, b"noodles", "noodles").await?;

        Ok(())
    }
}
