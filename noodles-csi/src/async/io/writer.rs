mod index;

use noodles_bgzf as bgzf;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use self::index::write_index;
use crate::Index;

/// An async CSI writer.
pub struct Writer<W> {
    inner: bgzf::r#async::io::Writer<W>,
}

impl<W> Writer<W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// use tokio::io;
    /// let writer = csi::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &bgzf::r#async::io::Writer<W> {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// use tokio::io;
    /// let mut writer = csi::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut bgzf::r#async::io::Writer<W> {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// use tokio::io;
    /// let writer = csi::r#async::io::Writer::new(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(self) -> bgzf::r#async::io::Writer<W> {
        self.inner
    }
}

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Creates an async CSI writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let writer = csi::r#async::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self {
            inner: bgzf::r#async::io::Writer::new(inner),
        }
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
    /// use noodles_csi as csi;
    /// let mut writer = csi::r#async::io::Writer::new(Vec::new());
    /// writer.shutdown().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn shutdown(&mut self) -> io::Result<()> {
        self.inner.shutdown().await
    }

    /// Writes a CSI index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_csi as csi;
    ///
    /// let index = csi::Index::default();
    ///
    /// let mut writer = csi::r#async::io::Writer::new(Vec::new());
    /// writer.write_index(&index).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_index(&mut self, index: &Index) -> io::Result<()> {
        write_index(&mut self.inner, index).await
    }
}
