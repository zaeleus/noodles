use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::fai::{Index, Record};

/// An async FASTA index (FAI) writer.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// use tokio::io;
    /// let writer = fai::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// use tokio::io;
    /// let mut writer = fai::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// use tokio::io;
    /// let writer = fai::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }
}

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Creates an async BAM index (BAI) writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::fai;
    /// use tokio::io;
    /// let writer = fai::r#async::io::Writer::new(io::sink());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Shuts down the output stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_fasta::fai;
    /// use tokio::io;
    /// let mut writer = fai::r#async::io::Writer::new(io::sink());
    /// writer.shutdown().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn shutdown(&mut self) -> io::Result<()> {
        self.inner.shutdown().await
    }

    /// Writes a FASTA index.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_fasta::fai;
    /// use tokio::io;
    ///
    /// let mut writer = fai::r#async::io::Writer::new(io::sink());
    ///
    /// let index = fai::Index::default();
    /// writer.write_index(&index).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_index(&mut self, index: &Index) -> io::Result<()> {
        for record in index.as_ref() {
            write_record(&mut self.inner, record).await?;
        }

        Ok(())
    }
}

async fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::fai::io::writer::write_record;

    let mut buf = Vec::new();
    write_record(&mut buf, record)?;
    writer.write_all(&buf).await
}
