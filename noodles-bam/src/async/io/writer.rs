mod header;

use noodles_bgzf as bgzf;
use noodles_sam as sam;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use self::header::write_header;
use crate::Record;

/// An async BAM writer.
pub struct Writer<W> {
    inner: W,
    buf: Vec<u8>,
}

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use tokio::io;
    /// let writer = bam::r#async::io::Writer::from(io::sink());
    /// let _inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use tokio::io;
    /// let mut writer = bam::r#async::io::Writer::from(io::sink());
    /// let _inner = writer.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use tokio::io;
    /// let writer = bam::r#async::io::Writer::from(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }

    /// Shuts down the output stream.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_bam as bam;
    /// use tokio::io;
    /// let mut writer = bam::r#async::io::Writer::new(io::sink());
    /// writer.shutdown().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn shutdown(&mut self) -> io::Result<()> {
        self.inner.shutdown().await
    }

    /// Writes a SAM header.
    ///
    /// This writes the BAM magic number, the raw SAM header, and a copy of the reference sequence
    /// dictionary as binary reference sequences.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_bam as bam;
    /// use noodles_sam as sam;
    /// use tokio::io;
    ///
    /// let mut writer = bam::r#async::io::Writer::new(io::sink());
    ///
    /// let header = sam::Header::default();
    /// writer.write_header(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        write_header(&mut self.inner, header).await
    }

    /// Writes a BAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_bam as bam;
    /// use noodles_sam as sam;
    /// use tokio::io;
    ///
    /// let mut writer = bam::r#async::io::Writer::new(io::sink());
    ///
    /// let header = sam::Header::default();
    /// let record = bam::Record::default();
    /// writer.write_record(&header, &record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_record(&mut self, header: &sam::Header, record: &Record) -> io::Result<()> {
        self.write_alignment_record(header, record).await
    }

    /// Writes an alignment record.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> tokio::io::Result<()> {
    /// use noodles_bam as bam;
    /// use noodles_sam::{self as sam, alignment::RecordBuf};
    /// use tokio::io;
    ///
    /// let mut writer = bam::r#async::io::Writer::new(io::sink());
    ///
    /// let header = sam::Header::default();
    /// let record = RecordBuf::default();
    /// writer.write_alignment_record(&header, &record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_alignment_record(
        &mut self,
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
    ) -> io::Result<()> {
        use crate::record::codec::encode;

        self.buf.clear();
        encode(&mut self.buf, header, record)?;

        let block_size = u32::try_from(self.buf.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        self.inner.write_u32_le(block_size).await?;

        self.inner.write_all(&self.buf).await?;

        Ok(())
    }
}

impl<W> Writer<bgzf::r#async::io::Writer<W>>
where
    W: AsyncWrite + Unpin,
{
    /// Creates an async BAM writer with a default compression level.
    ///
    /// The given stream is wrapped in a BGZF encoder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// use tokio::io;
    /// let writer = bam::r#async::io::Writer::new(io::sink());
    /// ```
    pub fn new(inner: W) -> Self {
        Self::from(bgzf::r#async::io::Writer::new(inner))
    }
}

impl<W> From<W> for Writer<W> {
    fn from(inner: W) -> Self {
        Self {
            inner,
            buf: Vec::new(),
        }
    }
}
