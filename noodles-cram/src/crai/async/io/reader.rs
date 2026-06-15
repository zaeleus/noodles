mod record;

use async_compression::tokio::bufread::GzipDecoder;
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, BufReader};

use self::record::read_record;
use crate::crai::{Index, Record};

/// An async CRAM index reader.
pub struct Reader<R> {
    inner: BufReader<GzipDecoder<BufReader<R>>>,
    buf: String,
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
            buf: String::new(),
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
        read_index(&mut self.inner, &mut self.buf).await
    }

    /// Reads a record.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_cram::crai;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.cram.crai")
    ///     .await
    ///     .map(crai::r#async::io::Reader::new)?;
    ///
    /// let mut record = crai::Record::default();
    ///
    /// while reader.read_record(&mut record).await? != 0 {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        self.buf.clear();
        read_record(&mut self.inner, &mut self.buf, record).await
    }
}

async fn read_index<R>(reader: &mut R, buf: &mut String) -> io::Result<Index>
where
    R: AsyncBufRead + Unpin,
{
    let mut index = Vec::new();
    let mut record = Record::default();

    while read_record(reader, buf, &mut record).await? != 0 {
        index.push(record.clone());
    }

    Ok(index)
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
    use noodles_core::Position;

    use super::*;
    use crate::crai::Record;

    #[tokio::test]
    async fn test_read_index() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"0\t10946\t6765\t17711\t233\t317811\n";

        let mut reader = &data[..];
        let mut buf = String::new();
        let actual = read_index(&mut reader, &mut buf).await?;

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

    #[tokio::test]
    async fn test_read_line() -> io::Result<()> {
        async fn t(buf: &mut String, mut src: &[u8], expected: &str) -> io::Result<()> {
            buf.clear();
            read_line(&mut src, buf).await?;
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
