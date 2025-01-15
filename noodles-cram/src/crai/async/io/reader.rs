use async_compression::tokio::bufread::GzipDecoder;
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, BufReader};

use crate::crai::Index;

/// An async CRAM index reader.
pub struct Reader<R> {
    inner: BufReader<GzipDecoder<BufReader<R>>>,
}

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    /// Creates an async CRAM index reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::crai;
    /// let data = [];
    /// let reader = crai::r#async::io::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self {
            inner: BufReader::new(GzipDecoder::new(BufReader::new(inner))),
        }
    }

    /// Reads a CRAM index.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_cram::crai;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.cram.crai")
    ///     .await
    ///     .map(crai::r#async::io::Reader::new)?;
    ///
    /// let index = reader.read_index().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_index(&mut self) -> io::Result<Index> {
        read_index(&mut self.inner).await
    }
}

async fn read_index<R>(reader: &mut R) -> io::Result<Index>
where
    R: AsyncBufRead + Unpin,
{
    let mut index = Vec::new();
    let mut buf = String::new();

    loop {
        buf.clear();

        match read_line(reader, &mut buf).await {
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

async fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    const LINE_FEED: char = '\n';

    match reader.read_line(buf).await {
        Ok(0) => Ok(0),
        Ok(n) => {
            if buf.ends_with(LINE_FEED) {
                buf.pop();
            }

            Ok(n)
        }
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;

    use super::*;
    use crate::crai::Record;

    #[tokio::test]
    async fn test_read_index() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"0\t10946\t6765\t17711\t233\t317811\n";

        let mut reader = &data[..];
        let actual = read_index(&mut reader).await?;

        let expected = vec![Record::new(
            Some(0),
            Position::new(10946),
            6765,
            17711,
            233,
            317811,
        )];

        assert_eq!(actual, expected);

        Ok(())
    }
}
