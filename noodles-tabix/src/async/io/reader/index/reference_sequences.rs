mod bins;
mod intervals;
mod metadata;

use noodles_csi::binning_index::index::{
    ReferenceSequence, reference_sequence::index::LinearIndex,
};
use tokio::io::{self, AsyncRead};

use self::{bins::read_bins, intervals::read_intervals, metadata::read_metadata};

pub(super) async fn read_reference_sequences<R>(
    reader: &mut R,
    reference_sequence_count: usize,
) -> io::Result<Vec<ReferenceSequence<LinearIndex>>>
where
    R: AsyncRead + Unpin,
{
    let mut reference_sequences = Vec::with_capacity(reference_sequence_count);

    for _ in 0..reference_sequence_count {
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
