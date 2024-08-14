//! Async alignment reader.
//!
//! Constructing a [`Reader`] is best done via a builder using the [`Builder`] type.

mod builder;

use futures::{Stream, StreamExt};
use noodles_bam as bam;
use noodles_cram as cram;
use noodles_sam as sam;
use std::pin::Pin;
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
    /// Reads the SAM header
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::alignment::r#async::io::reader::Builder;
    ///
    /// let data = b"@HD\tVN:1.6
    /// @SQ\tSN:chr1\tLN:2489
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// ";
    ///
    /// let mut reader = Builder::default().build_from_reader(&data[..]).await?;
    /// let header = reader.read_header().await?;
    ///
    /// assert_eq!(header.reference_sequences().len(), 1);
    /// # Ok(())
    /// # }
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
    /// use noodles_sam::alignment::Record;
    ///
    /// let data = b"@HD\tVN:1.6
    /// @SQ\tSN:chr1\tLN:2489
    /// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
    /// chr1\t0\tr1\t1\t60\t100M\t*\t0\t0\t*\t*
    /// ";
    ///
    /// let mut reader = Builder::default().build_from_reader(&data[..]).await?;
    /// let header = reader.read_header().await?;
    /// let mut records = reader.records(&header);
    ///
    /// let mut num_unmapped = 0;
    /// while let Some(record) = records.try_next().await? {
    ///     let is_unmapped = record.flags().map(|f| f.is_unmapped())?;
    ///     if is_unmapped {
    ///        num_unmapped += 1;
    ///    }
    /// }
    /// assert_eq!(num_unmapped, 1);
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
