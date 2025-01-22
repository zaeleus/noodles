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
