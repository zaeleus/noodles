use tokio::io::{
    self, AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncSeek, AsyncSeekExt, SeekFrom,
};

/// An async FASTA reader.
pub struct Reader<R> {
    inner: R,
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
    /// let mut reader = fasta::AsyncReader::new(&data[..]);
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
    /// let mut reader = fasta::AsyncReader::new(&data[..]);
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
    /// let mut reader = fasta::AsyncReader::new(&data[..]);
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

async fn read_sequence<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    use memchr::memchr;

    use crate::reader::DEFINITION_PREFIX;

    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    let mut n = 0;

    loop {
        let reader_buf = reader.fill_buf().await?;

        if reader_buf
            .first()
            .map(|&b| b == DEFINITION_PREFIX)
            .unwrap_or(true)
        {
            break;
        }

        let len = match memchr(LINE_FEED, reader_buf) {
            Some(i) => {
                let line = &reader_buf[..i];

                if line.ends_with(&[CARRIAGE_RETURN]) {
                    let end = line.len() - 1;
                    buf.extend_from_slice(&line[..end]);
                } else {
                    buf.extend_from_slice(line);
                }

                i + 1
            }
            None => {
                buf.extend(reader_buf);
                reader_buf.len()
            }
        };

        reader.consume(len);

        n += len;
    }

    Ok(n)
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
    async fn test_read_sequence() -> io::Result<()> {
        async fn t(buf: &mut Vec<u8>, mut reader: &[u8], expected: &[u8]) -> io::Result<()> {
            buf.clear();
            read_sequence(&mut reader, buf).await?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, b"ACGT\n", b"ACGT").await?;
        t(&mut buf, b"ACGT\n>sq1\n", b"ACGT").await?;
        t(&mut buf, b"NNNN\nNNNN\nNN\n", b"NNNNNNNNNN").await?;

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
