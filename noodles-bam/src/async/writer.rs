use std::ffi::CString;

use noodles_bgzf as bgzf;
use noodles_sam::{
    self as sam,
    alignment::Record,
    header::record::value::{map::ReferenceSequence, Map},
};
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::writer::record::encode_record;

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
    /// let writer = bam::AsyncWriter::from(Vec::new());
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
    /// use noodles_bam as bam;
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
    /// use noodles_bam as bam;
    /// let writer = bam::AsyncWriter::from(Vec::new());
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
    /// use noodles_bam as bam;
    /// let mut writer = bam::AsyncWriter::new(Vec::new());
    /// writer.shutdown().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn shutdown(&mut self) -> io::Result<()> {
        self.inner.shutdown().await
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
    /// use noodles_bam as bam;
    /// use noodles_sam as sam;
    ///
    /// let mut writer = bam::AsyncWriter::new(Vec::new());
    ///
    /// let header = sam::Header::builder().add_comment("noodles-bam").build();
    /// writer.write_header(&header).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        write_header(&mut self.inner, header).await
    }

    /// Writes the binary reference sequences after the SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use noodles_bam as bam;
    /// use noodles_sam::{
    ///     self as sam,
    ///     header::record::value::{map::ReferenceSequence, Map},
    /// };
    ///
    /// let mut writer = bam::AsyncWriter::new(Vec::new());
    ///
    /// let header = sam::Header::builder()
    ///     .add_reference_sequence(Map::<ReferenceSequence>::new("sq0".parse()?, 8)?)
    ///     .add_comment("noodles-bam")
    ///     .build();
    ///
    /// writer.write_header(&header).await?;
    /// writer.write_reference_sequences(header.reference_sequences()).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_reference_sequences(
        &mut self,
        reference_sequences: &sam::header::ReferenceSequences,
    ) -> io::Result<()> {
        write_reference_sequences(&mut self.inner, reference_sequences).await
    }

    /// Writes a BAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_bam as bam;
    /// use noodles_sam::{self as sam, alignment::Record};
    ///
    /// let mut writer = bam::AsyncWriter::new(Vec::new());
    ///
    /// let header = sam::Header::default();
    /// let record = Record::default();
    /// writer.write_record(&header, &record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_record(&mut self, header: &sam::Header, record: &Record) -> io::Result<()> {
        self.buf.clear();
        encode_record(&mut self.buf, header, record)?;

        let block_size = u32::try_from(self.buf.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        self.inner.write_u32_le(block_size).await?;

        self.inner.write_all(&self.buf).await?;

        Ok(())
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
    /// use noodles_bam as bam;
    /// use noodles_sam::{self as sam, alignment::Record};
    ///
    /// let mut writer = bam::AsyncWriter::new(Vec::new());
    ///
    /// let header = sam::Header::default();
    /// let record = Record::default();
    /// writer.write_alignment_record(&header, &record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_alignment_record(
        &mut self,
        header: &sam::Header,
        record: &Record,
    ) -> io::Result<()> {
        self.write_record(header, record).await
    }
}

impl<W> Writer<bgzf::AsyncWriter<W>>
where
    W: AsyncWrite + Unpin,
{
    /// Creates an async BAM writer with a default compression level.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let writer = bam::AsyncWriter::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self::from(bgzf::AsyncWriter::new(inner))
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

async fn write_header<W>(writer: &mut W, header: &sam::Header) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::MAGIC_NUMBER;

    writer.write_all(MAGIC_NUMBER).await?;

    let text = header.to_string();
    let l_text =
        u32::try_from(text.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(l_text).await?;

    writer.write_all(text.as_bytes()).await?;

    Ok(())
}

async fn write_reference_sequences<W>(
    writer: &mut W,
    reference_sequences: &sam::header::ReferenceSequences,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let n_ref = u32::try_from(reference_sequences.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(n_ref).await?;

    for reference_sequence in reference_sequences.values() {
        write_reference_sequence(writer, reference_sequence).await?;
    }

    Ok(())
}

async fn write_reference_sequence<W>(
    writer: &mut W,
    reference_sequence: &Map<ReferenceSequence>,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let c_name = CString::new(reference_sequence.name().as_bytes())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let name = c_name.as_bytes_with_nul();

    let l_name =
        u32::try_from(name.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(l_name).await?;
    writer.write_all(name).await?;

    let l_ref = u32::try_from(usize::from(reference_sequence.length()))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(l_ref).await?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_reference_sequence() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();
        let reference_sequence = Map::<ReferenceSequence>::new("sq0".parse()?, 8)?;
        write_reference_sequence(&mut buf, &reference_sequence).await?;

        let expected = [
            0x04, 0x00, 0x00, 0x00, // l_name = 4
            0x73, 0x71, 0x30, 0x00, // name = b"sq0\x00"
            0x08, 0x00, 0x00, 0x00, // l_ref = 8
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
