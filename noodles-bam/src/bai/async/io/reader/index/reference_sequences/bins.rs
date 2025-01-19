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
    use crate::bai::DEPTH;

    const METADATA_ID: usize = Bin::metadata_id(DEPTH);

    fn duplicate_bin_error(id: usize) -> io::Result<(IndexMap<usize, Bin>, Option<Metadata>)> {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("duplicate bin ID: {id}"),
        ))
    }

    let n_bin = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut bins = IndexMap::with_capacity(n_bin);
    let mut metadata = None;

    for _ in 0..n_bin {
        let id = reader.read_u32_le().await.and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        if id == METADATA_ID {
            let m = read_metadata(reader).await?;

            if metadata.replace(m).is_some() {
                return duplicate_bin_error(id);
            }
        } else {
            let chunks = read_chunks(reader).await?;
            let bin = Bin::new(chunks);

            if bins.insert(id, bin).is_some() {
                return duplicate_bin_error(id);
            }
        }
    }

    Ok((bins, metadata))
}
