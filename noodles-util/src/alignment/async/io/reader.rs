//! Async alignment reader.

mod builder;
mod inner;

use futures::Stream;
use noodles_sam as sam;
use tokio::io::{self, AsyncRead};

pub use self::builder::Builder;
use self::inner::Inner;

/// An async alignment reader.
pub struct Reader<R>(Inner<R>)
where
    R: AsyncRead + Unpin;

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    /// Reads the SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::alignment::r#async::io::reader::Builder;
    /// use tokio::io;
    /// let mut reader = Builder::default().build_from_reader(io::empty()).await?;
    /// let _header = reader.read_header().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_header(&mut self) -> io::Result<sam::Header> {
        self.0.read_header().await
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_util::alignment::r#async::io::reader::Builder;
    /// use tokio::io;
    ///
    /// let mut reader = Builder::default().build_from_reader(io::empty()).await?;
    /// let header = reader.read_header().await?;
    ///
    /// let mut records = reader.records(&header);
    ///
    /// while let Some(record) = records.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
    ) -> impl Stream<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'r {
        self.0.records(header)
    }
}
