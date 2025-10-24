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
    /// Creates an async alignment reader.
    ///
    /// This attempts to autodetect the compression method and format of the input.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::alignment;
    /// use tokio::io;
    /// let reader = alignment::r#async::io::Reader::new(io::empty()).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn new(reader: R) -> io::Result<Self> {
        Builder::default().build_from_reader(reader).await
    }

    /// Reads the SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::alignment;
    /// use tokio::io;
    /// let mut reader = alignment::r#async::io::Reader::new(io::empty()).await?;
    /// let header = reader.read_header().await?;
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
    /// use noodles_util::alignment;
    /// use tokio::io;
    ///
    /// let mut reader = alignment::r#async::io::Reader::new(io::empty()).await?;
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
