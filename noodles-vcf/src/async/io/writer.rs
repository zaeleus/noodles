use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::{variant::io::Write, Header, Record};

/// An async VCF writer.
///
/// If the inner writer is buffered, a call to [`Self::shutdown`] must be made before the writer is
/// dropped.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Creates an async VCF writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let writer = vcf::r#async::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let writer = vcf::r#async::io::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let mut writer = vcf::r#async::io::Writer::new(Vec::new());
    /// assert!(writer.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let writer = vcf::r#async::io::Writer::new(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }

    /// Shuts down the output stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_vcf as vcf;
    /// let mut writer = vcf::r#async::io::Writer::new(Vec::new());
    /// writer.shutdown().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn shutdown(&mut self) -> io::Result<()> {
        self.inner.shutdown().await
    }

    /// Writes a VCF header.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_vcf as vcf;
    ///
    /// let mut writer = vcf::r#async::io::Writer::new(Vec::new());
    ///
    /// let header = vcf::Header::default();
    /// writer.write_header(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_header(&mut self, header: &Header) -> io::Result<()> {
        let mut writer = crate::io::Writer::new(Vec::new());
        writer.write_header(header)?;
        self.inner.write_all(writer.get_ref()).await
    }

    /// Writes a VCF record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_core::Position;
    /// use noodles_vcf as vcf;
    ///
    /// let mut writer = vcf::r#async::io::Writer::new(Vec::new());
    ///
    /// let header = vcf::Header::default();
    /// let record = vcf::Record::default();
    /// writer.write_record(&header, &record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_record(&mut self, header: &Header, record: &Record) -> io::Result<()> {
        let mut writer = crate::io::Writer::new(Vec::new());
        writer.write_record(header, record)?;
        self.inner.write_all(writer.get_ref()).await?;
        Ok(())
    }

    /// Writes a variant record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_core::Position;
    /// use noodles_vcf as vcf;
    ///
    /// let header = vcf::Header::default();
    ///
    /// let record = vcf::variant::RecordBuf::builder()
    ///     .set_reference_sequence_name("sq0")
    ///     .set_variant_start(Position::MIN)
    ///     .set_reference_bases("A")
    ///     .build();
    ///
    /// let mut writer = vcf::r#async::io::Writer::new(Vec::new());
    /// writer.write_variant_record(&header, &record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_variant_record(
        &mut self,
        header: &Header,
        record: &dyn crate::variant::Record,
    ) -> io::Result<()> {
        let mut writer = crate::io::Writer::new(Vec::new());
        writer.write_variant_record(header, record)?;
        self.inner.write_all(writer.get_ref()).await?;
        Ok(())
    }
}
