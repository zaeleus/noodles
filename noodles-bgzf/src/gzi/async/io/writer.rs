mod index;

use tokio::io::{self, AsyncWrite};

use self::index::write_index;
use crate::gzi::Index;

/// An async GZ index (GZI) writer.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::gzi;
    /// use tokio::io;
    /// let writer = gzi::r#async::io::Writer::new(io::sink());
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
    /// use noodles_bgzf::gzi;
    /// use tokio::io;
    /// let mut writer = gzi::r#async::io::Writer::new(io::sink());
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
    /// use noodles_bgzf::gzi;
    /// use tokio::io;
    /// let writer = gzi::r#async::io::Writer::new(io::sink());
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
    /// Creates an async GZ index writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf::gzi;
    /// use tokio::io;
    /// let writer = gzi::r#async::io::Writer::new(io::sink());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Writes a GZ index.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_bgzf::gzi;
    /// use tokio::io;
    /// let index = gzi::Index::default();
    /// let mut writer = gzi::r#async::io::Writer::new(io::sink());
    /// writer.write_index(&index).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_index(&mut self, index: &Index) -> io::Result<()> {
        write_index(&mut self.inner, index).await
    }
}
