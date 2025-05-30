mod bins;
mod metadata;

use tokio::io::{self, AsyncRead, AsyncReadExt};

use self::{bins::read_bins, metadata::read_metadata};
use crate::binning_index::index::{ReferenceSequence, reference_sequence::index::BinnedIndex};

pub(super) async fn read_reference_sequences<R>(
    reader: &mut R,
    depth: u8,
) -> io::Result<Vec<ReferenceSequence<BinnedIndex>>>
where
    R: AsyncRead + Unpin,
{
    let n_ref = reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut reference_sequences = Vec::with_capacity(n_ref);

    for _ in 0..n_ref {
        let reference_sequence = read_reference_sequence(reader, depth).await?;
        reference_sequences.push(reference_sequence);
    }

    Ok(reference_sequences)
}

async fn read_reference_sequence<R>(
    reader: &mut R,
    depth: u8,
) -> io::Result<ReferenceSequence<BinnedIndex>>
where
    R: AsyncRead + Unpin,
{
    let (bins, index, metadata) = read_bins(reader, depth).await?;
    Ok(ReferenceSequence::new(bins, index, metadata))
}
