use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::{alignment::RecordBuf, Header};

/// An async SAM writer.
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
    /// Creates an async SAM writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let writer = sam::AsyncWriter::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let writer = sam::AsyncWriter::new(Vec::new());
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
    /// use noodles_sam as sam;
    /// let mut writer = sam::AsyncWriter::new(Vec::new());
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
    /// use noodles_sam as sam;
    /// let writer = sam::AsyncWriter::new(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }

    /// Writes a SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_sam as sam;
    ///
    /// let mut writer = sam::AsyncWriter::new(Vec::new());
    ///
    /// let header = sam::Header::builder().add_comment("noodles-sam").build();
    /// writer.write_header(&header).await?;
    ///
    /// assert_eq!(writer.get_ref(), b"@CO\tnoodles-sam\n");
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_header(&mut self, header: &Header) -> io::Result<()> {
        let mut writer = crate::io::Writer::new(Vec::new());
        writer.write_header(header)?;
        self.inner.write_all(writer.get_ref()).await
    }

    /// Writes a SAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_sam::{self as sam, alignment::RecordBuf};
    ///
    /// let mut writer = sam::AsyncWriter::new(Vec::new());
    ///
    /// let header = sam::Header::default();
    /// let record = RecordBuf::default();
    /// writer.write_record(&header, &record).await?;
    ///
    /// assert_eq!(writer.get_ref(), b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_record(&mut self, header: &Header, record: &RecordBuf) -> io::Result<()> {
        self.write_alignment_record(header, record).await
    }

    /// Writes an alignment record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_sam::{self as sam, alignment::RecordBuf};
    ///
    /// let mut writer = sam::AsyncWriter::new(Vec::new());
    ///
    /// let header = sam::Header::default();
    /// let record = RecordBuf::default();
    /// writer.write_alignment_record(&header, &record).await?;
    ///
    /// assert_eq!(writer.get_ref(), b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_alignment_record(
        &mut self,
        header: &Header,
        record: &RecordBuf,
    ) -> io::Result<()> {
        use crate::io::writer::write_record;

        let mut buf = Vec::new();
        write_record(&mut buf, header, record)?;
        self.inner.write_all(&buf).await
    }
}
