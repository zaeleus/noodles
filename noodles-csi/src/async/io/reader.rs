mod index;

use noodles_bgzf as bgzf;
use tokio::io::{self, AsyncRead};

use self::index::read_index;
use crate::Index;

/// An async CSI reader.
pub struct Reader<R>
where
    R: AsyncRead,
{
    inner: bgzf::AsyncReader<R>,
}

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    /// Creates an async CSI reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let data = [];
    /// let reader = csi::r#async::io::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self {
            inner: bgzf::AsyncReader::new(inner),
        }
    }

    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// use tokio::io;
    /// let reader = csi::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &bgzf::AsyncReader<R> {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// use tokio::io;
    /// let mut reader = csi::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut bgzf::AsyncReader<R> {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// use tokio::io;
    /// let reader = csi::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> bgzf::AsyncReader<R> {
        self.inner
    }

    /// Reads the CSI index.
    ///
    /// The position of the stream is expected to be at the beginning.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_csi as csi;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bcf.csi")
    ///     .await
    ///     .map(csi::r#async::io::Reader::new)?;
    ///
    /// let index = reader.read_index().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_index(&mut self) -> io::Result<Index> {
        read_index(&mut self.inner).await
    }
}
