//! Async alignment writer.

mod builder;
mod inner;

use noodles_sam as sam;
use tokio::io::{self, AsyncWrite};

pub use self::builder::Builder;
use self::inner::Inner;

/// An async alignment writer.
pub struct Writer<W>(Inner<W>)
where
    W: AsyncWrite;

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Writes a SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_sam as sam;
    /// use noodles_util::alignment::r#async::io::writer::Builder;
    /// use tokio::io::{self, AsyncWriteExt};
    ///
    /// let mut writer = Builder::default().build_from_writer(io::sink()).await?;
    ///
    /// let header = sam::Header::default();
    /// writer.write_header(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        self.0.write_header(header).await
    }

    /// Writes an alignment record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_sam as sam;
    /// use noodles_util::alignment::r#async::io::writer::Builder;
    /// use tokio::io::{self, AsyncWriteExt};
    ///
    /// let mut writer = Builder::default().build_from_writer(io::sink()).await?;
    ///
    /// let header = sam::Header::default();
    /// let record = sam::Record::default();
    /// writer.write_record(&header, &record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_record(
        &mut self,
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
    ) -> io::Result<()> {
        self.0.write_record(header, record).await
    }

    /// Shuts down the alignment format writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_sam as sam;
    /// use noodles_util::alignment::{self, io::Format};
    /// use tokio::io;
    ///
    /// let mut writer = alignment::r#async::io::writer::Builder::default()
    ///     .set_format(Format::Sam)
    ///     .build_from_writer(io::sink()).await?;
    ///
    /// let header = sam::Header::default();
    /// writer.shutdown(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn shutdown(&mut self, header: &sam::Header) -> io::Result<()> {
        self.0.shutdown(header).await
    }
}
