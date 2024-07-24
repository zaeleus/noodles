use noodles_bcf as bcf;
use noodles_vcf as vcf;
use tokio::io::{self, AsyncWrite};

/// An async variant writer.
pub enum Writer<W> {
    /// BCF.
    Bcf(bcf::r#async::io::Writer<W>),
    /// VCF.
    Vcf(vcf::r#async::io::Writer<W>),
}

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
    /// use noodles_util::variant::r#async::io::Writer;
    /// use noodles_vcf as vcf;
    /// use tokio::io;
    ///
    /// let mut writer = Writer::Vcf(vcf::r#async::io::Writer::new(io::sink()));
    ///
    /// let header = vcf::Header::default();
    /// writer.write_header(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_header(&mut self, header: &vcf::Header) -> io::Result<()> {
        match self {
            Self::Bcf(writer) => writer.write_header(header).await,
            Self::Vcf(writer) => writer.write_header(header).await,
        }
    }

    /// Writes a variant record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// #  async fn main() -> tokio::io::Result<()> {
    /// use noodles_util::variant::r#async::io::Writer;
    /// use noodles_vcf as vcf;
    /// use tokio::io;
    ///
    /// let mut writer = Writer::Vcf(vcf::r#async::io::Writer::new(io::sink()));
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
        match self {
            Self::Bcf(writer) => writer.write_variant_record(header, record).await,
            Self::Vcf(writer) => writer.write_variant_record(header, record).await,
        }
    }
}
