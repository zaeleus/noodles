use async_compression::tokio::write::GzipEncoder;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::crai::Record;

/// An async CRAM index writer.
pub struct Writer<W> {
    inner: GzipEncoder<W>,
}

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Creates an async CRAM index writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::crai;
    /// let writer = crai::r#async::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self {
            inner: GzipEncoder::new(inner),
        }
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::crai;
    /// let writer = crai::r#async::io::Writer::new(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> W {
        self.inner.into_inner()
    }

    /// Shuts down the output stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_cram::crai;
    /// let mut writer = crai::r#async::io::Writer::new(Vec::new());
    /// writer.shutdown().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn shutdown(&mut self) -> io::Result<()> {
        self.inner.shutdown().await
    }

    /// Writes a CRAM index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_core::Position;
    /// use noodles_cram::crai;
    ///
    /// let mut writer = crai::r#async::io::Writer::new(Vec::new());
    ///
    /// let index = vec![crai::Record::new(
    ///     Some(0),
    ///     Position::new(10946),
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// )];
    ///
    /// writer.write_index(&index).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_index(&mut self, index: &[Record]) -> io::Result<()> {
        write_index(&mut self.inner, index).await
    }
}

async fn write_index<W>(writer: &mut W, index: &[Record]) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    for record in index {
        write_record(writer, record).await?;
    }

    Ok(())
}

async fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    const LINE_FEED: u8 = b'\n';

    writer.write_all(record.to_string().as_bytes()).await?;
    writer.write_all(&[LINE_FEED]).await?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_record() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_core::Position;

        let mut buf = Vec::new();

        let record = Record::new(Some(0), Position::new(10946), 6765, 17711, 233, 317811);
        write_record(&mut buf, &record).await?;

        let expected = b"0\t10946\t6765\t17711\t233\t317811\n";
        assert_eq!(buf, expected);

        Ok(())
    }
}
