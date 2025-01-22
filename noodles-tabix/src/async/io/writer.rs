mod index;

use noodles_bgzf as bgzf;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use self::index::write_index;
use crate::Index;

/// An async tabix writer.
pub struct Writer<W> {
    inner: bgzf::AsyncWriter<W>,
}

impl<W> Writer<W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// use tokio::io;
    /// let writer = tabix::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &bgzf::AsyncWriter<W> {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// use tokio::io;
    /// let mut writer = tabix::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut bgzf::AsyncWriter<W> {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// use tokio::io;
    /// let writer = tabix::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(self) -> bgzf::AsyncWriter<W> {
        self.inner
    }
}

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Creates an async tabix writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let writer = tabix::r#async::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self {
            inner: bgzf::AsyncWriter::new(inner),
        }
    }

    /// Shuts down the output stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> std::io::Result<()> {
    /// use noodles_tabix as tabix;
    /// let mut writer = tabix::r#async::io::Writer::new(Vec::new());
    /// writer.shutdown().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn shutdown(&mut self) -> io::Result<()> {
        self.inner.shutdown().await
    }

    /// Writes a tabix index.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> std::io::Result<()> {
    /// use noodles_csi::binning_index::index::Header;
    /// use noodles_tabix as tabix;
    ///
    /// let index = tabix::Index::builder().set_header(Header::default()).build();
    ///
    /// let mut writer = tabix::r#async::io::Writer::new(Vec::new());
    /// writer.write_index(&index).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_index(&mut self, index: &Index) -> io::Result<()> {
        write_index(&mut self.inner, index).await
    }
}
