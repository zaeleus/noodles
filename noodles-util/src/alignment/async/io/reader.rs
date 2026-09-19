//! Async alignment reader.

mod builder;
mod inner;

use futures::Stream;
use noodles_core::Region;
use noodles_sam as sam;
use tokio::io::{self, AsyncRead, AsyncSeek};

pub use self::builder::Builder;
use self::inner::Inner;
use crate::alignment::Index;

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

impl<R> Reader<R>
where
    R: AsyncRead + AsyncSeek + Unpin,
{
    /// Returns a stream over records that intersects the given region.
    ///
    /// To query for unmapped records, use [`Self::query_unmapped`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use futures::TryStreamExt;
    /// use noodles_util::alignment;
    ///
    /// let mut reader = alignment::r#async::io::reader::Builder::default()
    ///     .build_from_path("sample.bam")
    ///     .await?;
    ///
    /// let header = reader.read_header().await?;
    ///
    /// let index = alignment::r#async::fs::read_associated_index("sample.bam").await?;
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(&header, &index, &region)?;
    ///
    /// while let Some(record) = query.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query<'r, 'h: 'r, 'i: 'r>(
        &'r mut self,
        header: &'h sam::Header,
        index: &'i Index,
        region: &Region,
    ) -> io::Result<impl Stream<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'r> {
        self.0.query(header, index, region)
    }

    /// Returns a stream of unmapped records after querying for the unmapped region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use futures::TryStreamExt;
    /// use noodles_util::alignment;
    ///
    /// let mut reader = alignment::r#async::io::reader::Builder::default()
    ///     .build_from_path("sample.bam")
    ///     .await?;
    ///
    /// let header = reader.read_header().await?;
    ///
    /// let index = alignment::r#async::fs::read_associated_index("sample.bam").await?;
    /// let mut query = reader.query_unmapped(&header, &index).await?;
    ///
    /// while let Some(record) = query.try_next().await? {
    ///     // ...
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub async fn query_unmapped<'r, 'h: 'r, 'i: 'r>(
        &'r mut self,
        header: &'h sam::Header,
        index: &'i Index,
    ) -> io::Result<impl Stream<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'r> {
        self.0.query_unmapped(header, index).await
    }
}
