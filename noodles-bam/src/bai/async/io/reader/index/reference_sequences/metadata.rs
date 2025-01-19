use noodles_bgzf as bgzf;
use noodles_csi::binning_index::index::reference_sequence::Metadata;
use tokio::io::{self, AsyncRead, AsyncReadExt};

pub(super) async fn read_metadata<R>(reader: &mut R) -> io::Result<Metadata>
where
    R: AsyncRead + Unpin,
{
    const METADATA_CHUNK_COUNT: usize = 2;

    let n_chunk = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    if n_chunk != METADATA_CHUNK_COUNT {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "invalid metadata pseudo-bin chunk count: expected {METADATA_CHUNK_COUNT}, got {n_chunk}"
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

#[cfg(test)]
mod tests {
    use super::*;

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
}
