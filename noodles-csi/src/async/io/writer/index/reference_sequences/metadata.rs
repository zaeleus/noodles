use noodles_bgzf as bgzf;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use crate::binning_index::index::reference_sequence::{Bin, Metadata};

pub(super) async fn write_metadata<W>(
    writer: &mut W,
    depth: u8,
    metadata: &Metadata,
) -> io::Result<()>
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
