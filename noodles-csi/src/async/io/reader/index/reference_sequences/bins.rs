mod chunks;

use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use tokio::io::{self, AsyncRead, AsyncReadExt};

use self::chunks::read_chunks;
use super::read_metadata;
use crate::binning_index::index::reference_sequence::{Bin, Metadata, index::BinnedIndex};

pub(super) async fn read_bins<R>(
    reader: &mut R,
    depth: u8,
) -> io::Result<(IndexMap<usize, Bin>, BinnedIndex, Option<Metadata>)>
where
    R: AsyncRead + Unpin,
{
    fn duplicate_bin_error(
        id: usize,
    ) -> io::Result<(IndexMap<usize, Bin>, BinnedIndex, Option<Metadata>)> {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("duplicate bin ID: {id}"),
        ))
    }

    let n_bin = reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut bins = IndexMap::with_capacity(n_bin);
    let mut index = BinnedIndex::with_capacity(n_bin);

    let metadata_id = Bin::metadata_id(depth);
    let mut metadata = None;

    for _ in 0..n_bin {
        let id = reader.read_u32_le().await.and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let loffset = reader
            .read_u64_le()
            .await
            .map(bgzf::VirtualPosition::from)?;

        if id == metadata_id {
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

            index.insert(id, loffset);
        }
    }

    Ok((bins, index, metadata))
}
