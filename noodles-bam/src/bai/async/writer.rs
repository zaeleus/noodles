use noodles_bgzf as bgzf;
use noodles_csi::{
    binning_index::ReferenceSequenceExt,
    index::reference_sequence::{bin::Chunk, Metadata},
    BinningIndex,
};
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::bai::{
    index::reference_sequence::{Bin, ReferenceSequence},
    Index, MAGIC_NUMBER,
};

/// An async BAM index (BAI) writer.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Creates an async BAM index (BAI) writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// let writer = bai::AsyncWriter::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// let writer = bai::AsyncWriter::new(Vec::new());
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
    /// use noodles_bam::bai;
    /// let mut writer = bai::AsyncWriter::new(Vec::new());
    /// writer.shutdown().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn shutdown(&mut self) -> io::Result<()> {
        self.inner.shutdown().await
    }

    /// Writes a BAM index (BAI) header.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_bam::bai;
    /// let mut writer = bai::AsyncWriter::new(Vec::new());
    /// writer.write_header().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_header(&mut self) -> io::Result<()> {
        write_magic(&mut self.inner).await
    }

    /// Writes a BAM index.
    ///
    /// The position of the stream is expected to be directly after the header.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_bam::bai;
    ///
    /// let index = bai::Index::default();
    ///
    /// let mut writer = bai::AsyncWriter::new(Vec::new());
    /// writer.write_header().await?;
    /// writer.write_index(&index).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_index(&mut self, index: &Index) -> io::Result<()> {
        write_reference_sequences(&mut self.inner, index.reference_sequences()).await?;

        if let Some(unplaced_unmapped_record_count) = index.unplaced_unmapped_record_count() {
            write_unplaced_unmapped_record_count(&mut self.inner, unplaced_unmapped_record_count)
                .await?;
        }

        Ok(())
    }
}

async fn write_magic<W>(writer: &mut W) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    writer.write_all(MAGIC_NUMBER).await
}

async fn write_reference_sequences<W>(
    writer: &mut W,
    reference_sequences: &[ReferenceSequence],
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let n_ref = u32::try_from(reference_sequences.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(n_ref).await?;

    for reference_sequence in reference_sequences {
        write_reference_sequence(writer, reference_sequence).await?;
    }

    Ok(())
}

async fn write_reference_sequence<W>(
    writer: &mut W,
    reference_sequence: &ReferenceSequence,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    write_bins(
        writer,
        reference_sequence.bins(),
        reference_sequence.metadata(),
    )
    .await?;

    write_intervals(writer, reference_sequence.intervals()).await?;

    Ok(())
}

async fn write_bins<W>(writer: &mut W, bins: &[Bin], metadata: Option<&Metadata>) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let n_bin = u32::try_from(bins.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        .and_then(|n| {
            if metadata.is_some() {
                n.checked_add(1)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "n_bin overflow"))
            } else {
                Ok(n)
            }
        })?;

    writer.write_u32_le(n_bin).await?;

    for bin in bins {
        write_bin(writer, bin).await?;
    }

    if let Some(m) = metadata {
        write_metadata(writer, m).await?;
    }

    Ok(())
}

async fn write_bin<W>(writer: &mut W, bin: &Bin) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let id = u32::try_from(bin.id()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(id).await?;
    write_chunks(writer, bin.chunks()).await?;
    Ok(())
}

async fn write_chunks<W>(writer: &mut W, chunks: &[Chunk]) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let n_chunk =
        u32::try_from(chunks.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(n_chunk).await?;

    for chunk in chunks {
        write_chunk(writer, chunk).await?;
    }

    Ok(())
}

async fn write_chunk<W>(writer: &mut W, chunk: &Chunk) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let chunk_beg = u64::from(chunk.start());
    writer.write_u64_le(chunk_beg).await?;

    let chunk_end = u64::from(chunk.end());
    writer.write_u64_le(chunk_end).await?;

    Ok(())
}

async fn write_intervals<W>(writer: &mut W, intervals: &[bgzf::VirtualPosition]) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let n_intv = u32::try_from(intervals.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(n_intv).await?;

    for &interval in intervals {
        let ioffset = u64::from(interval);
        writer.write_u64_le(ioffset).await?;
    }

    Ok(())
}

async fn write_metadata<W>(writer: &mut W, metadata: &Metadata) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::bai::index::reference_sequence::bin::{METADATA_CHUNK_COUNT, METADATA_ID};

    let id =
        u32::try_from(METADATA_ID).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(id).await?;

    let n_chunk = u32::try_from(METADATA_CHUNK_COUNT)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(n_chunk).await?;

    let ref_beg = u64::from(metadata.start_position());
    writer.write_u64_le(ref_beg).await?;

    let ref_end = u64::from(metadata.end_position());
    writer.write_u64_le(ref_end).await?;

    let n_mapped = metadata.mapped_record_count();
    writer.write_u64_le(n_mapped).await?;

    let n_unmapped = metadata.unmapped_record_count();
    writer.write_u64_le(n_unmapped).await?;

    Ok(())
}

async fn write_unplaced_unmapped_record_count<W>(
    writer: &mut W,
    unplaced_unmapped_record_count: u64,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    writer.write_u64_le(unplaced_unmapped_record_count).await
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_magic() -> io::Result<()> {
        let mut buf = Vec::new();
        write_magic(&mut buf).await?;
        assert_eq!(buf, b"BAI\x01");
        Ok(())
    }

    #[tokio::test]
    async fn test_write_bins() -> io::Result<()> {
        let bins = vec![Bin::new(8, Vec::new())];

        let mut buf = Vec::new();
        write_bins(&mut buf, &bins, None).await?;

        let expected = [
            0x01, 0x00, 0x00, 0x00, // n_bins = 1
            0x08, 0x00, 0x00, 0x00, // bins[0].bin = 8
            0x00, 0x00, 0x00, 0x00, // bins[0].n_chunk = 0
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[tokio::test]
    async fn test_write_bins_with_metadata() -> io::Result<()> {
        let bins = vec![Bin::new(8, Vec::new())];
        let metadata = Metadata::new(
            bgzf::VirtualPosition::from(13),
            bgzf::VirtualPosition::from(21),
            5,
            0,
        );

        let mut buf = Vec::new();
        write_bins(&mut buf, &bins, Some(&metadata)).await?;

        #[rustfmt::skip]
        let expected = [
            0x02, 0x00, 0x00, 0x00, // n_bins = 2

            0x08, 0x00, 0x00, 0x00, // bins[0].bin = 8
            0x00, 0x00, 0x00, 0x00, // bins[0].n_chunk = 0

            0x4a, 0x92, 0x00, 0x00, // bins[1].bin = 37450
            0x02, 0x00, 0x00, 0x00, // bins[1].n_chunk = 2
            0x0d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].chunks[0].chunk_beg = 13
            0x15, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].chunks[0].chunk_end = 21
            0x05, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].chunks[1].chunk_beg = 5
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // bins[1].chunks[1].chunk_end = 0
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[tokio::test]
    async fn test_write_bin() -> io::Result<()> {
        let bin = Bin::new(
            8,
            vec![Chunk::new(
                bgzf::VirtualPosition::from(13),
                bgzf::VirtualPosition::from(21),
            )],
        );

        let mut buf = Vec::new();
        write_bin(&mut buf, &bin).await?;

        let expected = [
            0x08, 0x00, 0x00, 0x00, // bin = 8
            0x01, 0x00, 0x00, 0x00, // n_chunk = 1
            0x0d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // chunk[0].chunk_beg
            0x15, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // chunk[0].chunk_end
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[tokio::test]
    async fn test_write_metadata() -> io::Result<()> {
        let metadata = Metadata::new(
            bgzf::VirtualPosition::from(610),
            bgzf::VirtualPosition::from(1597),
            55,
            0,
        );

        let mut buf = Vec::new();
        write_metadata(&mut buf, &metadata).await?;

        let expected = [
            0x4a, 0x92, 0x00, 0x00, // bin = 37450
            0x02, 0x00, 0x00, 0x00, // n_chunks = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_unmapped = 0
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
