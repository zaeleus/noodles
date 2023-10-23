use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::{Header, Record};

const LINE_FEED: u8 = b'\n';

/// An async VCF writer.
pub struct Writer<W>
where
    W: AsyncWrite,
{
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
    /// let writer = vcf::AsyncWriter::new(Vec::new());
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
    /// let mut writer = vcf::AsyncWriter::new(Vec::new());
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
    /// let mut writer = bam::AsyncWriter::from(Vec::new());
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
    /// let writer = bam::AsyncWriter::from(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
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
    /// let mut writer = vcf::AsyncWriter::new(Vec::new());
    ///
    /// let header = vcf::Header::default();
    /// writer.write_header(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_header(&mut self, header: &Header) -> io::Result<()> {
        let raw_header = header.to_string();
        self.inner.write_all(raw_header.as_bytes()).await
    }

    /// Writes a VCF record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use noodles_vcf::{self as vcf, record::Position};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// let mut writer = vcf::AsyncWriter::new(Vec::new());
    /// writer.write_record(&record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_record(&mut self, record: &Record) -> io::Result<()> {
        let raw_record = record.to_string();
        self.inner.write_all(raw_record.as_bytes()).await?;
        self.inner.write_u8(LINE_FEED).await?;
        Ok(())
    }
}
