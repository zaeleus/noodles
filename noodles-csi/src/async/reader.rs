use noodles_bgzf as bgzf;
use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::{
    index::{
        reference_sequence::{bin::Chunk, Bin, Metadata},
        ReferenceSequence,
    },
    Index,
};

/// An async CSI reader.
pub struct Reader<R>
where
    R: AsyncRead,
{
    inner: bgzf::AsyncReader<R>,
}

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    /// Creates an async CSI reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_csi as csi;
    /// let data = [];
    /// let reader = csi::AsyncReader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self {
            inner: bgzf::AsyncReader::new(inner),
        }
    }

    /// Reads the CSI index.
    ///
    /// The position of the stream is expected to be at the beginning.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// #
    /// # #[tokio::main]
    /// # async fn main() -> io::Result<()> {
    /// use noodles_csi as csi;
    /// use tokio::fs::File;
    ///
    /// let mut reader = File::open("sample.bcf.csi")
    ///     .await
    ///     .map(csi::AsyncReader::new)?;
    ///
    /// let index = reader.read_index().await?;
    /// # Ok(())
    /// # }
    /// ```
    pub async fn read_index(&mut self) -> io::Result<Index> {
        read_magic(&mut self.inner).await?;

        let (min_shift, depth, aux) = read_header(&mut self.inner).await?;
        let reference_sequences = read_reference_sequences(&mut self.inner, depth).await?;
        let unplaced_unmapped_record_count =
            read_unplaced_unmapped_record_count(&mut self.inner).await?;

        let mut builder = Index::builder()
            .set_min_shift(min_shift)
            .set_depth(depth)
            .set_aux(aux)
            .set_reference_sequences(reference_sequences);

        if let Some(count) = unplaced_unmapped_record_count {
            builder = builder.set_unplaced_unmapped_record_count(count);
        }

        Ok(builder.build())
    }
}

async fn read_magic<R>(reader: &mut R) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    use crate::MAGIC_NUMBER;

    let mut magic = [0; 4];
    reader.read_exact(&mut magic).await?;

    if magic == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid CSI header",
        ))
    }
}

async fn read_header<R>(reader: &mut R) -> io::Result<(u8, u8, Vec<u8>)>
where
    R: AsyncRead + Unpin,
{
    let min_shift = reader
        .read_i32_le()
        .await
        .and_then(|n| u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))?;

    let depth = reader
        .read_i32_le()
        .await
        .and_then(|n| u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))?;

    let aux = read_aux(reader).await?;

    Ok((min_shift, depth, aux))
}

async fn read_aux<R>(reader: &mut R) -> io::Result<Vec<u8>>
where
    R: AsyncRead + Unpin,
{
    let l_aux = reader.read_i32_le().await.and_then(|len| {
        usize::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut aux = vec![0; l_aux];
    reader.read_exact(&mut aux).await?;

    Ok(aux)
}

async fn read_reference_sequences<R>(
    reader: &mut R,
    depth: u8,
) -> io::Result<Vec<ReferenceSequence>>
where
    R: AsyncRead + Unpin,
{
    let n_ref = reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut reference_sequences = Vec::with_capacity(n_ref);

    for _ in 0..n_ref {
        let reference_sequence = read_reference_sequence(reader, depth).await?;
        reference_sequences.push(reference_sequence);
    }

    Ok(reference_sequences)
}

async fn read_reference_sequence<R>(reader: &mut R, depth: u8) -> io::Result<ReferenceSequence>
where
    R: AsyncRead + Unpin,
{
    let (bins, metadata) = read_bins(reader, depth).await?;
    Ok(ReferenceSequence::new(bins, metadata))
}

async fn read_bins<R>(reader: &mut R, depth: u8) -> io::Result<(Vec<Bin>, Option<Metadata>)>
where
    R: AsyncRead + Unpin,
{
    let n_bin = reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut bins = Vec::with_capacity(n_bin);

    let metadata_id = Bin::metadata_id(depth);
    let mut metadata = None;

    for _ in 0..n_bin {
        let id = reader.read_u32_le().await.and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let loffset = reader
            .read_u64_le()
            .await
            .map(bgzf::VirtualPosition::from)?;

        if id == metadata_id {
            metadata = read_metadata(reader).await.map(Some)?;
        } else {
            let chunks = read_chunks(reader).await?;
            let bin = Bin::new(id, loffset, chunks);
            bins.push(bin);
        }
    }

    Ok((bins, metadata))
}

async fn read_chunks<R>(reader: &mut R) -> io::Result<Vec<Chunk>>
where
    R: AsyncRead + Unpin,
{
    let n_chunk = reader.read_i32_le().await.and_then(|n| {
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

async fn read_metadata<R>(reader: &mut R) -> io::Result<Metadata>
where
    R: AsyncRead + Unpin,
{
    use crate::index::reference_sequence::bin::METADATA_CHUNK_COUNT;

    let n_chunk = reader.read_u32_le().await?;

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
        let data = b"CSI\x01";
        let mut reader = &data[..];
        assert!(read_magic(&mut reader).await.is_ok());

        let data = [];
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
    async fn test_read_header() -> io::Result<()> {
        let data = [
            0x0e, 0x00, 0x00, 0x00, // min_shift = 14
            0x05, 0x00, 0x00, 0x00, // depth = 5
            0x04, 0x00, 0x00, 0x00, // l_aux = 4
            0x6e, 0x64, 0x6c, 0x73, // aux = b"ndls"
        ];

        let mut reader = &data[..];
        let (min_shift, depth, aux) = read_header(&mut reader).await?;

        assert_eq!(min_shift, 14);
        assert_eq!(depth, 5);
        assert_eq!(aux, b"ndls");

        Ok(())
    }
}
