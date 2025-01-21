mod chunks;

use indexmap::IndexMap;
use noodles_csi::binning_index::index::reference_sequence::{Bin, Metadata};
use tokio::io::{self, AsyncRead, AsyncReadExt};

use self::chunks::read_chunks;
use super::read_metadata;

pub(super) async fn read_bins<R>(
    reader: &mut R,
) -> io::Result<(IndexMap<usize, Bin>, Option<Metadata>)>
where
    R: AsyncRead + Unpin,
{
    use crate::index::DEPTH;

    const METADATA_ID: usize = Bin::metadata_id(DEPTH);

    let n_bin = reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut bins = IndexMap::with_capacity(n_bin);
    let mut metadata = None;

    for _ in 0..n_bin {
        let id = reader.read_u32_le().await.and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let is_duplicate = if id == METADATA_ID {
            let m = read_metadata(reader).await?;
            metadata.replace(m).is_some()
        } else {
            let chunks = read_chunks(reader).await?;
            let bin = Bin::new(chunks);
            bins.insert(id, bin).is_some()
        };

        if is_duplicate {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("duplicate bin ID: {id}"),
            ));
        }
    }

    Ok((bins, metadata))
}
