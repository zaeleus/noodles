//! Async alignment writer.

mod builder;

use noodles_bam as bam;
use noodles_cram as cram;
use noodles_sam as sam;
use tokio::io::{self, AsyncWrite};

pub use self::builder::Builder;

/// An async alignment writer.
pub enum Writer<W: tokio::io::AsyncWrite> {
    /// SAM.
    Sam(sam::r#async::io::Writer<W>),
    /// BAM.
    Bam(bam::r#async::io::Writer<W>),
    /// CRAM.
    Cram(cram::r#async::io::Writer<W>),
}

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
    /// use tokio::io;
    /// use tokio::io::AsyncWriteExt;
    ///
    /// let mut writer = Builder::default().build_from_writer(io::sink());
    ///
    /// let header = sam::Header::default();
    /// writer.write_header(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        match self {
            Self::Sam(writer) => writer.write_header(header).await,
            Self::Bam(writer) => writer.write_header(header).await,
            Self::Cram(writer) => writer.write_file_header(header).await,
        }
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
    /// use tokio::io;
    /// use tokio::io::AsyncWriteExt;
    ///
    /// let mut writer = Builder::default().build_from_writer(io::sink());
    ///
    /// let record = sam::Record::default();
    /// writer.write_record(&record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_record(
        &mut self,
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
    ) -> io::Result<()> {
        match self {
            Self::Sam(writer) => writer.write_alignment_record(header, record).await,
            Self::Bam(writer) => writer.write_alignment_record(header, record).await,
            Self::Cram(writer) => writer.write_record(header, record).await,
        }
    }
}
