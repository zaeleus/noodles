//! Async alignment reader.

mod builder;

use std::pin::Pin;

use futures::{Stream, StreamExt};
use noodles_bam as bam;
use noodles_cram as cram;
use noodles_sam as sam;
use tokio::io::{self, AsyncBufRead};

pub use self::builder::Builder;

/// An async alignment reader.
pub enum Reader<R> {
    /// SAM.
    Sam(sam::r#async::io::Reader<R>),
    /// BAM.
    Bam(bam::r#async::io::Reader<R>),
    /// CRAM.
    Cram(cram::r#async::io::Reader<R>),
}

impl<R> Reader<R>
where
    R: AsyncBufRead + Unpin,
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
        match self {
            Self::Sam(reader) => reader.read_header().await,
            Self::Bam(reader) => reader.read_header().await,
            Self::Cram(reader) => reader.read_header().await,
        }
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
        #[allow(clippy::type_complexity)]
        let records: Pin<
            Box<dyn Stream<Item = io::Result<Box<dyn sam::alignment::Record>>>>,
        > = match self {
            Self::Sam(reader) => Box::pin(
                reader
                    .records()
                    .map(|result| result.map(|r| Box::new(r) as Box<dyn sam::alignment::Record>)),
            ),
            Self::Bam(reader) => Box::pin(
                reader
                    .records()
                    .map(|result| result.map(|r| Box::new(r) as Box<dyn sam::alignment::Record>)),
            ),
            Self::Cram(reader) => Box::pin(
                reader
                    .records(header)
                    .map(|result| result.map(|r| Box::new(r) as Box<dyn sam::alignment::Record>)),
            ),
        };

        records
    }
}
