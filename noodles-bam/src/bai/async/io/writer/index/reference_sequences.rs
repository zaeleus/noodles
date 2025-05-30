mod bins;
mod intervals;
mod metadata;

use noodles_csi::binning_index::{
    ReferenceSequence as _,
    index::{ReferenceSequence, reference_sequence::index::LinearIndex},
};
use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use self::{bins::write_bins, intervals::write_intervals, metadata::write_metadata};

pub(super) async fn write_reference_sequences<W>(
    writer: &mut W,
    reference_sequences: &[ReferenceSequence<LinearIndex>],
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    let n_ref = u32::try_from(reference_sequences.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_u32_le(n_ref).await?;

    for reference_sequence in reference_sequences {
        write_reference_sequence(writer, reference_sequence).await?;
    }

    Ok(())
}

async fn write_reference_sequence<W>(
    writer: &mut W,
    reference_sequence: &ReferenceSequence<LinearIndex>,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    write_bins(
        writer,
        reference_sequence.bins(),
        reference_sequence.metadata(),
    )
    .await?;

    write_intervals(writer, reference_sequence.index()).await?;

    Ok(())
}
