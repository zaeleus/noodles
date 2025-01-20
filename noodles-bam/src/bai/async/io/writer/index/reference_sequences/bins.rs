mod chunks;

use indexmap::IndexMap;
use noodles_csi::binning_index::index::reference_sequence::{Bin, Metadata};
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use self::chunks::write_chunks;
use super::write_metadata;

pub(super) async fn write_bins<W>(
    writer: &mut W,
    bins: &IndexMap<usize, Bin>,
    metadata: Option<&Metadata>,
) -> io::Result<()>
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

    for (&id, bin) in bins {
        write_bin(writer, id, bin).await?;
    }

    if let Some(m) = metadata {
        write_metadata(writer, m).await?;
    }

    Ok(())
}

async fn write_bin<W>(writer: &mut W, id: usize, bin: &Bin) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let id = u32::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(id).await?;
    write_chunks(writer, bin.chunks()).await?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_bgzf as bgzf;
    use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;

    use super::*;

    #[tokio::test]
    async fn test_write_bins() -> io::Result<()> {
        let bins = [(8, Bin::new(Vec::new()))].into_iter().collect();

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
        let bins = [(8, Bin::new(Vec::new()))].into_iter().collect();
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
        let bin = Bin::new(vec![Chunk::new(
            bgzf::VirtualPosition::from(13),
            bgzf::VirtualPosition::from(21),
        )]);

        let mut buf = Vec::new();
        write_bin(&mut buf, 8, &bin).await?;

        let expected = [
            0x08, 0x00, 0x00, 0x00, // bin = 8
            0x01, 0x00, 0x00, 0x00, // n_chunk = 1
            0x0d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // chunk[0].chunk_beg
            0x15, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // chunk[0].chunk_end
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
