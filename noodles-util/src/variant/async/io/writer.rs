//! Async variant writer.

mod builder;
mod inner;

use noodles_vcf as vcf;
use tokio::io::{self, AsyncWrite};

pub use self::builder::Builder;
use self::inner::Inner;

/// An async variant writer.
pub struct Writer<W>(Inner<W>)
where
    W: AsyncWrite;

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Writes a VCF header.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// #  async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::variant::r#async::io::writer::Builder;
    /// use noodles_vcf as vcf;
    /// use tokio::io;
    ///
    /// let mut writer = Builder::default().build_from_writer(io::sink());
    ///
    /// let header = vcf::Header::default();
    /// writer.write_header(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_header(&mut self, header: &vcf::Header) -> io::Result<()> {
        self.0.write_header(header).await
    }

    /// Writes a variant record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// #  async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::variant::r#async::io::writer::Builder;
    /// use noodles_vcf as vcf;
    /// use tokio::io;
    ///
    /// let mut writer = Builder::default().build_from_writer(io::sink());
    ///
    /// let header = vcf::Header::default();
    /// writer.write_header(&header).await?;
    ///
    /// let record = vcf::Record::default();
    /// writer.write_record(&header, &record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_record(
        &mut self,
        header: &vcf::Header,
        record: &dyn vcf::variant::Record,
    ) -> io::Result<()> {
        self.0.write_record(header, record).await
    }
}
