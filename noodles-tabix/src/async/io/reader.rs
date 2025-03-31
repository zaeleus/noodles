mod index;

use noodles_bgzf as bgzf;
use tokio::io::{self, AsyncRead};

use self::index::read_index;
use crate::Index;

/// An async tabix reader.
pub struct Reader<R>
where
    R: AsyncRead,
{
    inner: bgzf::r#async::io::Reader<R>,
}

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    /// Creates an async tabix reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// let data = [];
    /// let reader = tabix::r#async::io::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self {
            inner: bgzf::r#async::io::Reader::new(inner),
        }
    }

    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// use tokio::io;
    /// let reader = tabix::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &bgzf::r#async::io::Reader<R> {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// use tokio::io;
    /// let mut reader = tabix::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut bgzf::r#async::io::Reader<R> {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_tabix as tabix;
    /// use tokio::io;
    /// let reader = tabix::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> bgzf::r#async::io::Reader<R> {
        self.inner
    }

    /// Reads the tabix index.
    ///
    /// The position of the stream is expected to be at the beginning.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_tabix as tabix;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.vcf.gz.tbi")
    ///     .await
    ///     .map(tabix::r#async::io::Reader::new)?;
    ///
    /// let index = reader.read_index().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_index(&mut self) -> io::Result<Index> {
        read_index(&mut self.inner).await
    }
}
