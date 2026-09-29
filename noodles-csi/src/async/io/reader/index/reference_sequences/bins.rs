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
    let bin_count = read_bin_count(reader).await?;

    let mut bins = IndexMap::with_capacity(bin_count);
    let mut index = BinnedIndex::with_capacity(bin_count);

    let metadata_id = Bin::metadata_id(depth)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid depth"))?;

    let mut metadata = None;

    for _ in 0..bin_count {
        let id = reader.read_u32_le().await.and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let loffset = reader
            .read_u64_le()
            .await
            .map(bgzf::VirtualPosition::from)?;

        let is_duplicate = if id == metadata_id {
            let m = read_metadata(reader).await?;
            metadata.replace(m).is_some()
        } else {
            index.insert(id, loffset);

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

    Ok((bins, index, metadata))
}

async fn read_bin_count<R>(reader: &mut R) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    reader
        .read_i32_le()
        .await
        .and_then(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_bin_count() -> io::Result<()> {
        assert_eq!(read_bin_count(&mut &[0x00, 0x00, 0x00, 0x00][..]).await?, 0);
        assert_eq!(read_bin_count(&mut &[0x08, 0x00, 0x00, 0x00][..]).await?, 8);

        #[cfg(not(target_pointer_width = "16"))]
        assert_eq!(
            read_bin_count(&mut &[0xff, 0xff, 0xff, 0x7f][..]).await?,
            i32::MAX as usize
        );

        assert!(matches!(
            read_bin_count(&mut io::empty()).await,
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        assert!(matches!(
            read_bin_count(&mut &[0xff, 0xff, 0xff, 0xff][..]).await,
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
