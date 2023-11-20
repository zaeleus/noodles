use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::{
    binning_index::ReferenceSequence as _,
    index::{
        reference_sequence::{bin::Chunk, index::BinnedIndex, Bin, Metadata},
        Header, ReferenceSequence,
    },
    BinningIndex, Index,
};

/// An async CSI writer.
pub struct Writer<W> {
    inner: bgzf::AsyncWriter<W>,
}

impl<W> Writer<W>
where
    W: AsyncWrite + Unpin,
{
    /// Creates an async CSI writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let writer = csi::AsyncWriter::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self {
            inner: bgzf::AsyncWriter::new(inner),
        }
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let writer = csi::AsyncWriter::new(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> W {
        self.inner.into_inner()
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
    /// use noodles_csi as csi;
    /// let mut writer = csi::AsyncWriter::new(Vec::new());
    /// writer.shutdown().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn shutdown(&mut self) -> io::Result<()> {
        self.inner.shutdown().await
    }

    /// Writes a CSI index.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_csi as csi;
    ///
    /// let index = csi::Index::default();
    ///
    /// let mut writer = csi::AsyncWriter::new(Vec::new());
    /// writer.write_index(&index).await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn write_index(&mut self, index: &Index<BinnedIndex>) -> io::Result<()> {
        write_index(&mut self.inner, index).await
    }
}

async fn write_index<W>(writer: &mut W, index: &Index<BinnedIndex>) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    write_magic(writer).await?;

    write_header(writer, index.min_shift(), index.depth(), index.header()).await?;
    write_reference_sequences(writer, index.depth(), index.reference_sequences()).await?;

    if let Some(unplaced_unmapped_record_count) = index.unplaced_unmapped_record_count() {
        writer.write_u64_le(unplaced_unmapped_record_count).await?;
    }

    Ok(())
}

async fn write_magic<W>(writer: &mut W) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::io::MAGIC_NUMBER;

    writer.write_all(MAGIC_NUMBER).await
}

async fn write_header<W>(
    writer: &mut W,
    min_shift: u8,
    depth: u8,
    header: Option<&Header>,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    writer.write_i32_le(i32::from(min_shift)).await?;
    writer.write_i32_le(i32::from(depth)).await?;
    write_aux(writer, header).await?;
    Ok(())
}

async fn write_aux<W>(writer: &mut W, header: Option<&Header>) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::writer::index::header::write_header as write_tabix_header;

    let mut aux = Vec::new();

    if let Some(hdr) = header {
        write_tabix_header(&mut aux, hdr)?;
    }

    let l_aux =
        i32::try_from(aux.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32_le(l_aux).await?;

    writer.write_all(&aux).await?;

    Ok(())
}

async fn write_reference_sequences<W>(
    writer: &mut W,
    depth: u8,
    reference_sequences: &[ReferenceSequence<BinnedIndex>],
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    for reference_sequence in reference_sequences {
        write_reference_sequence(writer, depth, reference_sequence).await?;
    }

    Ok(())
}

async fn write_reference_sequence<W>(
    writer: &mut W,
    depth: u8,
    reference_sequence: &ReferenceSequence<BinnedIndex>,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    write_bins(
        writer,
        depth,
        reference_sequence.bins(),
        reference_sequence.index(),
        reference_sequence.metadata(),
    )
    .await
}

async fn write_bins<W>(
    writer: &mut W,
    depth: u8,
    bins: &IndexMap<usize, Bin>,
    index: &IndexMap<usize, bgzf::VirtualPosition>,
    metadata: Option<&Metadata>,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let n_bin = i32::try_from(bins.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        .and_then(|n| {
            if metadata.is_some() {
                n.checked_add(1)
                    .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "n_bin overflow"))
            } else {
                Ok(n)
            }
        })?;

    writer.write_i32_le(n_bin).await?;

    for (id, bin) in bins {
        let first_record_start_position = index.get(id).copied().unwrap_or_default();
        write_bin(writer, *id, first_record_start_position, bin).await?;
    }

    if let Some(m) = metadata {
        write_metadata(writer, depth, m).await?;
    }

    Ok(())
}

async fn write_bin<W>(
    writer: &mut W,
    id: usize,
    first_record_start_position: bgzf::VirtualPosition,
    bin: &Bin,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let bin_id = u32::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(bin_id).await?;

    let loffset = u64::from(first_record_start_position);
    writer.write_u64_le(loffset).await?;

    write_chunks(writer, bin.chunks()).await?;

    Ok(())
}

async fn write_chunks<W>(writer: &mut W, chunks: &[Chunk]) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let n_chunk =
        i32::try_from(chunks.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32_le(n_chunk).await?;

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

async fn write_metadata<W>(writer: &mut W, depth: u8, metadata: &Metadata) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    const N_CHUNK: i32 = 2;

    let bin_id = u32::try_from(Bin::metadata_id(depth))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(bin_id).await?;

    let loffset = u64::from(bgzf::VirtualPosition::default());
    writer.write_u64_le(loffset).await?;

    writer.write_i32_le(N_CHUNK).await?;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_magic() -> io::Result<()> {
        let mut buf = Vec::new();
        write_magic(&mut buf).await?;
        assert_eq!(buf, b"CSI\x01");
        Ok(())
    }

    #[tokio::test]
    async fn test_write_header() -> io::Result<()> {
        let mut buf = Vec::new();
        write_header(&mut buf, 14, 5, None).await?;

        let expected = [
            0x0e, 0x00, 0x00, 0x00, // min_shift = 14
            0x05, 0x00, 0x00, 0x00, // depth = 5
            0x00, 0x00, 0x00, 0x00, // l_aux = 0
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[tokio::test]
    async fn test_write_metadata() -> io::Result<()> {
        let mut buf = Vec::new();
        let depth = 5;
        let metadata = Metadata::new(
            bgzf::VirtualPosition::from(610),
            bgzf::VirtualPosition::from(1597),
            55,
            0,
        );

        write_metadata(&mut buf, depth, &metadata).await?;

        let expected = [
            0x4a, 0x92, 0x00, 0x00, // bin = 37450
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // loffset = 0
            0x02, 0x00, 0x00, 0x00, // chunks = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_unmapped = 0
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
