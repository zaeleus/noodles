use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::binning_index::index::reference_sequence::bin::Chunk;

pub(super) async fn write_chunks<W>(writer: &mut W, chunks: &[Chunk]) -> io::Result<()>
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
