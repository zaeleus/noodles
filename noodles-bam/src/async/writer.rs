mod builder;
mod record;

pub use self::builder::Builder;

use std::{ffi::CString, mem};

use noodles_bgzf as bgzf;
use noodles_sam as sam;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use self::record::write_sam_record;
use crate::Record;

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
    /// Creates an async BAM writer builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    /// let builder = bam::AsyncWriter::builder(Vec::new());
    /// let writer = builder.build();
    /// ```
    pub fn builder(inner: W) -> Builder<W> {
        Builder::new(inner)
    }

    /// Creates an async BAM writer with a default compression level.
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
    /// let mut writer = bam::AsyncWriter::new(Vec::new());
    /// let record = bam::Record::default();
    /// writer.write_record(&record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_record(&mut self, record: &Record) -> io::Result<()> {
        write_record(&mut self.inner, record).await
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
    /// use noodles_bam as bam;
    /// use noodles_sam as sam;
    ///
    /// let mut writer = bam::AsyncWriter::new(Vec::new());
    ///
    /// let reference_sequences = sam::header::ReferenceSequences::default();
    /// let record = sam::Record::default();
    /// writer.write_sam_record(&reference_sequences, &record).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_sam_record(
        &mut self,
        reference_sequences: &sam::header::ReferenceSequences,
        record: &sam::Record,
    ) -> io::Result<()> {
        write_sam_record(&mut self.inner, reference_sequences, record).await
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

async fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let block_size = u32::try_from(record.block_size())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(block_size).await?;

    writer.write_i32_le(record.ref_id).await?;
    writer.write_i32_le(record.pos).await?;

    let l_read_name = u8::try_from(record.read_name.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u8(l_read_name).await?;

    let mapq = u8::from(record.mapping_quality());
    writer.write_u8(mapq).await?;

    writer.write_u16_le(record.bin).await?;

    let n_cigar_op = u16::try_from(record.cigar().len() / mem::size_of::<u32>())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u16_le(n_cigar_op).await?;

    let flag = u16::from(record.flags());
    writer.write_u16_le(flag).await?;

    let l_seq = u32::try_from(record.sequence().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(l_seq).await?;

    writer.write_i32_le(record.next_ref_id).await?;
    writer.write_i32_le(record.next_pos).await?;

    writer.write_i32_le(record.tlen).await?;

    writer.write_all(&record.read_name).await?;

    for &raw_op in record.cigar().iter() {
        writer.write_u32_le(raw_op).await?;
    }

    writer.write_all(record.sequence().as_ref()).await?;
    writer.write_all(record.quality_scores()).await?;
    writer.write_all(record.data()).await?;

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
