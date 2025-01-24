mod bins;
mod metadata;

use tokio::io::{self, AsyncWrite};

use self::{bins::write_bins, metadata::write_metadata};
use crate::binning_index::{
    index::{reference_sequence::index::BinnedIndex, ReferenceSequence},
    ReferenceSequence as _,
};

pub(super) async fn write_reference_sequences<W>(
    writer: &mut W,
    depth: u8,
    reference_sequences: &[ReferenceSequence<BinnedIndex>],
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    for reference_sequence in reference_sequences {
        write_reference_sequence(writer, depth, reference_sequence).await?;
    }

    Ok(())
}

async fn write_reference_sequence<W>(
    writer: &mut W,
    depth: u8,
    reference_sequence: &ReferenceSequence<BinnedIndex>,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    write_bins(
        writer,
        depth,
        reference_sequence.bins(),
        reference_sequence.index(),
        reference_sequence.metadata(),
    )
    .await
}
