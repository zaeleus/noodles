mod index;

use tokio::io::{self, AsyncRead};

use self::index::read_index;
use crate::bai::Index;

/// An async BAM index (BAI) reader.
pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R> {
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// use tokio::io;
    /// let reader = bai::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.get_ref();
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// use tokio::io;
    /// let mut reader = bai::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// use tokio::io;
    /// let reader = bai::r#async::io::Reader::new(io::empty());
    /// let _inner = reader.into_inner();
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    /// Creates an async BAM index (BAI) reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// use tokio::io;
    /// let reader = bai::r#async::io::Reader::new(io::empty());
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads the BAM index.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_bam::bai;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam.bai")
    ///     .await
    ///     .map(bai::r#async::io::Reader::new)?;
    ///
    /// let index = reader.read_index().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_index(&mut self) -> io::Result<Index> {
        read_index(&mut self.inner).await
    }
}
