use noodles_bgzf as bgzf;
use noodles_csi::index::reference_sequence::{bin::Chunk, Metadata};
use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::bai::{
    self,
    index::{reference_sequence::Bin, ReferenceSequence},
    Index, MAGIC_NUMBER,
};

/// An async BAM index (BAI) reader.
pub struct Reader<R>
where
    R: AsyncRead,
{
    inner: R,
}

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    /// Creates an async BAM index (BAI) reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::bai;
    /// let data = [];
    /// let reader = bai::AsyncReader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads the BAM index (BAI) header.
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
    ///
    /// let data = b"BAI\x01";
    /// let mut reader = bai::AsyncReader::new(&data[..]);
    /// reader.read_header().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_header(&mut self) -> io::Result<()> {
        read_magic(&mut self.inner).await
    }

    /// Reads the BAM index.
    ///
    /// The position of the stream is expected to be directly after the header.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_bam::bai;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bam.bai")
    ///     .await
    ///     .map(bai::AsyncReader::new)?;
    ///
    /// reader.read_header().await?;
    ///
    /// let index = reader.read_index().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_index(&mut self) -> io::Result<Index> {
        let reference_sequences = read_reference_sequences(&mut self.inner).await?;

        let unplaced_unmapped_record_count =
            read_unplaced_unmapped_record_count(&mut self.inner).await?;

        Ok(Index::new(
            reference_sequences,
            unplaced_unmapped_record_count,
        ))
    }
}

async fn read_magic<R>(reader: &mut R) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    let mut magic = [0; 4];
    reader.read_exact(&mut magic).await?;

    if magic == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BAI header",
        ))
    }
}

async fn read_reference_sequences<R>(reader: &mut R) -> io::Result<Vec<ReferenceSequence>>
where
    R: AsyncRead + Unpin,
{
    let n_ref = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut reference_sequences = Vec::with_capacity(n_ref);

    for _ in 0..n_ref {
        let reference_sequence = read_reference_sequence(reader).await?;
        reference_sequences.push(reference_sequence);
    }

    Ok(reference_sequences)
}

async fn read_reference_sequence<R>(reader: &mut R) -> io::Result<ReferenceSequence>
where
    R: AsyncRead + Unpin,
{
    let (bins, metadata) = read_bins(reader).await?;
    let intervals = read_intervals(reader).await?;
    Ok(ReferenceSequence::new(bins, intervals, metadata))
}

async fn read_bins<R>(reader: &mut R) -> io::Result<(Vec<Bin>, Option<Metadata>)>
where
    R: AsyncRead + Unpin,
{
    use bai::index::reference_sequence::bin::METADATA_ID;

    let n_bin = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut bins = Vec::with_capacity(n_bin);
    let mut metadata = None;

    for _ in 0..n_bin {
        let id = reader.read_u32_le().await.and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        if id == METADATA_ID {
            metadata = read_metadata(reader).await.map(Some)?;
        } else {
            let chunks = read_chunks(reader).await?;
            let bin = Bin::new(id, chunks);
            bins.push(bin);
        }
    }

    Ok((bins, metadata))
}

async fn read_chunks<R>(reader: &mut R) -> io::Result<Vec<Chunk>>
where
    R: AsyncRead + Unpin,
{
    let n_chunk = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut chunks = Vec::with_capacity(n_chunk);

    for _ in 0..n_chunk {
        let chunk = read_chunk(reader).await?;
        chunks.push(chunk);
    }

    Ok(chunks)
}

async fn read_chunk<R>(reader: &mut R) -> io::Result<Chunk>
where
    R: AsyncRead + Unpin,
{
    let chunk_beg = reader
        .read_u64_le()
        .await
        .map(bgzf::VirtualPosition::from)?;

    let chunk_end = reader
        .read_u64_le()
        .await
        .map(bgzf::VirtualPosition::from)?;

    Ok(Chunk::new(chunk_beg, chunk_end))
}

async fn read_intervals<R>(reader: &mut R) -> io::Result<Vec<bgzf::VirtualPosition>>
where
    R: AsyncRead + Unpin,
{
    let n_intv = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut intervals = Vec::with_capacity(n_intv);

    for _ in 0..n_intv {
        let ioffset = reader
            .read_u64_le()
            .await
            .map(bgzf::VirtualPosition::from)?;

        intervals.push(ioffset);
    }

    Ok(intervals)
}

async fn read_metadata<R>(reader: &mut R) -> io::Result<Metadata>
where
    R: AsyncRead + Unpin,
{
    use bai::index::reference_sequence::bin::METADATA_CHUNK_COUNT;

    let n_chunk = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if n_chunk != METADATA_CHUNK_COUNT {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid metadata pseudo-bin chunk count: expected {}, got {}",
                METADATA_CHUNK_COUNT, n_chunk
            ),
        ));
    }

    let ref_beg = reader
        .read_u64_le()
        .await
        .map(bgzf::VirtualPosition::from)?;

    let ref_end = reader
        .read_u64_le()
        .await
        .map(bgzf::VirtualPosition::from)?;

    let n_mapped = reader.read_u64_le().await?;
    let n_unmapped = reader.read_u64_le().await?;

    Ok(Metadata::new(ref_beg, ref_end, n_mapped, n_unmapped))
}

async fn read_unplaced_unmapped_record_count<R>(reader: &mut R) -> io::Result<Option<u64>>
where
    R: AsyncRead + Unpin,
{
    match reader.read_u64_le().await {
        Ok(n_no_coor) => Ok(Some(n_no_coor)),
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(None),
        Err(e) => Err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_magic() {
        let data = b"BAI\x01";
        let mut reader = &data[..];
        assert!(read_magic(&mut reader).await.is_ok());
    }

    #[tokio::test]
    async fn test_read_magic_with_invalid_magic_number() {
        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"BAI";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"MThd";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader).await,
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[tokio::test]
    async fn test_read_chunk() -> io::Result<()> {
        let data = [
            0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // chunk_beg = 8
            0x0d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // chunk_end = 13
        ];

        let mut reader = &data[..];
        let actual = read_chunk(&mut reader).await?;

        let expected = Chunk::new(
            bgzf::VirtualPosition::from(8),
            bgzf::VirtualPosition::from(13),
        );

        assert_eq!(actual, expected);

        Ok(())
    }

    #[tokio::test]
    async fn test_read_intervals() -> io::Result<()> {
        let data = [
            0x03, 0x00, 0x00, 0x00, // n_intv = 3
            0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ioffset[0] = 8
            0x0d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ioffset[1] = 13
            0x15, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ioffset[2] = 21
        ];

        let mut reader = &data[..];
        let actual = read_intervals(&mut reader).await?;

        let expected = vec![
            bgzf::VirtualPosition::from(8),
            bgzf::VirtualPosition::from(13),
            bgzf::VirtualPosition::from(21),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[tokio::test]
    async fn test_read_metadata() -> io::Result<()> {
        let data = [
            0x02, 0x00, 0x00, 0x00, // n_chunk = 2
            0x62, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_beg = 610
            0x3d, 0x06, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // ref_end = 1597
            0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_mapped = 55
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // n_unmapped = 0
        ];

        let mut reader = &data[..];
        let actual = read_metadata(&mut reader).await?;

        let expected = Metadata::new(
            bgzf::VirtualPosition::from(610),
            bgzf::VirtualPosition::from(1597),
            55,
            0,
        );

        assert_eq!(actual, expected);

        Ok(())
    }

    #[tokio::test]
    async fn test_read_unplaced_unmapped_record_count() -> io::Result<()> {
        let data = [];
        let mut reader = &data[..];
        assert_eq!(
            read_unplaced_unmapped_record_count(&mut reader).await?,
            None
        );

        let data = [0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert_eq!(
            read_unplaced_unmapped_record_count(&mut reader).await?,
            Some(8)
        );

        Ok(())
    }
}
