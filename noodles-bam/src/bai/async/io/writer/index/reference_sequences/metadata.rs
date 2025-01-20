use noodles_csi::binning_index::index::reference_sequence::{Bin, Metadata};
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

pub(super) async fn write_metadata<W>(writer: &mut W, metadata: &Metadata) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    use crate::bai::DEPTH;

    const METADATA_ID: usize = Bin::metadata_id(DEPTH);
    const METADATA_CHUNK_COUNT: usize = 2;

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

#[cfg(test)]
mod tests {
    use noodles_bgzf as bgzf;

    use super::*;

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
