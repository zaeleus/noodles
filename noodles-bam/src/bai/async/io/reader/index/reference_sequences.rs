mod bins;
mod intervals;
mod metadata;

use noodles_csi::binning_index::index::{
    ReferenceSequence, reference_sequence::index::LinearIndex,
};
use tokio::io::{self, AsyncRead, AsyncReadExt};

use self::{bins::read_bins, intervals::read_intervals, metadata::read_metadata};

pub(super) async fn read_reference_sequences<R>(
    reader: &mut R,
) -> io::Result<Vec<ReferenceSequence<LinearIndex>>>
where
    R: AsyncRead + Unpin,
{
    let n_ref = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut reference_sequences = Vec::with_capacity(n_ref);

    for _ in 0..n_ref {
        let reference_sequence = read_reference_sequence(reader).await?;
        reference_sequences.push(reference_sequence);
    }

    Ok(reference_sequences)
}

async fn read_reference_sequence<R>(reader: &mut R) -> io::Result<ReferenceSequence<LinearIndex>>
where
    R: AsyncRead + Unpin,
{
    let (bins, metadata) = read_bins(reader).await?;
    let intervals = read_intervals(reader).await?;
    Ok(ReferenceSequence::new(bins, intervals, metadata))
}
