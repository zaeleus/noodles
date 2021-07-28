use std::{convert::TryFrom, ffi::CString};

use noodles_bgzf as bgzf;
use noodles_sam as sam;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::{Record, MAGIC_NUMBER};

/// An async BAM writer.
pub struct Writer<W>
where
    W: AsyncWrite,
{
    inner: bgzf::AsyncWriter<W>,
}

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Creates an async BAm writer with a default compression level.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let writer = bam::AsyncWriter::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self {
            inner: bgzf::AsyncWriter::new(inner),
        }
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
        self.inner.write_all(MAGIC_NUMBER).await?;

        let text = header.to_string();
        let l_text = u32::try_from(text.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        self.inner.write_u32_le(l_text).await?;

        self.inner.write_all(text.as_bytes()).await?;

        Ok(())
    }

    /// Writes the binary reference sequences after the SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// # #[tokio::main]
    /// # async fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// use noodles_bam as bam;
    /// use noodles_sam::{self as sam, header::ReferenceSequence};
    ///
    /// let mut writer = bam::AsyncWriter::new(Vec::new());
    ///
    /// let header = sam::Header::builder()
    ///     .add_reference_sequence(ReferenceSequence::new("sq0", 8)?)
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
        let n_ref = u32::try_from(reference_sequences.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        self.inner.write_u32_le(n_ref).await?;

        for reference_sequence in reference_sequences.values() {
            write_reference_sequence(&mut self.inner, reference_sequence).await?;
        }

        Ok(())
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
    /// let mut writer = bam::AsyncWriter::new(Vec::new());
    /// let record = bam::Record::default();
    /// writer.write_record(&record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_record(&mut self, record: &Record) -> io::Result<()> {
        let block_size = u32::try_from(record.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        self.inner.write_u32_le(block_size).await?;

        self.inner.write_all(record).await?;

        Ok(())
    }
}

async fn write_reference_sequence<W>(
    writer: &mut W,
    reference_sequence: &sam::header::ReferenceSequence,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let c_name = CString::new(reference_sequence.name())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let name = c_name.as_bytes_with_nul();

    let l_name =
        u32::try_from(name.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(l_name).await?;
    writer.write_all(name).await?;

    let l_ref = u32::try_from(reference_sequence.len())
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
        let reference_sequence = sam::header::ReferenceSequence::new("sq0", 8)?;
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
