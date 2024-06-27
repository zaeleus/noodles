mod sequence;

use tokio::io::{
    self, AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncSeek, AsyncSeekExt, SeekFrom,
};

use self::sequence::read_sequence;

/// An async FASTA reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// use tokio::io;
    /// let reader = fasta::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// use tokio::io;
    /// let mut reader = fasta::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// use tokio::io;
    /// let reader = fasta::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> Reader<R>
where
    R: AsyncBufRead + Unpin,
{
    /// Creates an async FASTA reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// let data = [];
    /// let mut reader = fasta::r#async::io::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads a raw definition line.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_fasta as fasta;
    ///
    /// let data = b">sq0\nACGT\n>sq1\nNNNN\nNNNN\nNN\n";
    /// let mut reader = fasta::r#async::io::Reader::new(&data[..]);
    ///
    /// let mut buf = String::new();
    /// reader.read_definition(&mut buf).await?;
    ///
    /// assert_eq!(buf, ">sq0");
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_definition(&mut self, buf: &mut String) -> io::Result<usize> {
        read_line(&mut self.inner, buf).await
    }

    /// Reads a sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_fasta as fasta;
    ///
    /// let data = b">sq0\nACGT\n>sq1\nNNNN\nNNNN\nNN\n";
    /// let mut reader = fasta::r#async::io::Reader::new(&data[..]);
    /// reader.read_definition(&mut String::new()).await?;
    ///
    /// let mut buf = Vec::new();
    /// reader.read_sequence(&mut buf).await?;
    ///
    /// assert_eq!(buf, b"ACGT");
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_sequence(&mut self, buf: &mut Vec<u8>) -> io::Result<usize> {
        read_sequence(&mut self.inner, buf).await
    }
}

impl<R> Reader<R>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    /// Seeks the underlying stream to the given position.
    pub async fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.inner.seek(pos).await
    }
}

pub(crate) async fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
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
    async fn test_read_definition() -> io::Result<()> {
        let data = b">sq0\nACGT\n";
        let mut reader = Reader::new(&data[..]);

        let mut buf = String::new();
        reader.read_definition(&mut buf).await?;

        assert_eq!(buf, ">sq0");

        Ok(())
    }

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
