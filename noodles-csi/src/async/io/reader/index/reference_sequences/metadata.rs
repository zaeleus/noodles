use noodles_bgzf as bgzf;
use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::binning_index::index::reference_sequence::{Metadata, bin::METADATA_CHUNK_COUNT};

pub(super) async fn read_metadata<R>(reader: &mut R) -> io::Result<Metadata>
where
    R: AsyncRead + Unpin,
{
    read_chunk_count(reader).await?;

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

async fn read_chunk_count<R>(reader: &mut R) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    let n = reader.read_u32_le().await?;

    if n == METADATA_CHUNK_COUNT {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid chunk count: expected {METADATA_CHUNK_COUNT}, got {n}"),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_chunk_count() {
        assert!(
            read_chunk_count(&mut &[0x02, 0x00, 0x00, 0x00][..])
                .await
                .is_ok()
        );

        assert!(matches!(
            read_chunk_count(&mut io::empty()).await,
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));
        assert!(matches!(
            read_chunk_count(&mut &[0x00, 0x00, 0x00, 0x00][..]).await,
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
        assert!(matches!(
            read_chunk_count(&mut &[0x08, 0x00, 0x00, 0x00][..]).await,
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }
}
