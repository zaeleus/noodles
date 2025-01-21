use noodles_bgzf as bgzf;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use tokio::io::{self, AsyncRead, AsyncReadExt};

pub(super) async fn read_chunks<R>(reader: &mut R) -> io::Result<Vec<Chunk>>
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
    let cnk_beg = reader
        .read_u64_le()
        .await
        .map(bgzf::VirtualPosition::from)?;

    let cnk_end = reader
        .read_u64_le()
        .await
        .map(bgzf::VirtualPosition::from)?;

    Ok(Chunk::new(cnk_beg, cnk_end))
}
