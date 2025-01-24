mod chunks;

use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use self::chunks::write_chunks;
use super::write_metadata;
use crate::binning_index::index::reference_sequence::{Bin, Metadata};

pub(super) async fn write_bins<W>(
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
