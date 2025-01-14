mod index;

use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use self::index::write_index;
use crate::bai::Index;

/// An async BAM index (BAI) writer.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// use tokio::io;
    /// let writer = bai::r#async::io::Writer::new(io::sink());
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
    /// use noodles_bam::bai;
    /// use tokio::io;
    /// let mut writer = bai::r#async::io::Writer::new(io::sink());
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
    /// use noodles_bam::bai;
    /// use tokio::io;
    /// let writer = bai::r#async::io::Writer::new(io::sink());
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
    /// use noodles_bam::bai;
    /// let writer = bai::r#async::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
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
    /// use noodles_bam::bai;
    /// let mut writer = bai::r#async::io::Writer::new(Vec::new());
    /// writer.shutdown().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn shutdown(&mut self) -> io::Result<()> {
        self.inner.shutdown().await
    }

    /// Writes a BAM index.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_bam::bai;
    ///
    /// let index = bai::Index::default();
    ///
    /// let mut writer = bai::r#async::io::Writer::new(Vec::new());
    /// writer.write_index(&index).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_index(&mut self, index: &Index) -> io::Result<()> {
        write_index(&mut self.inner, index).await
    }
}
