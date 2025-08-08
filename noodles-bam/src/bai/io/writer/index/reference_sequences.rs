mod bins;
mod intervals;
mod metadata;

use std::io::{self, Write};

use noodles_csi::binning_index::{
    ReferenceSequence as _,
    index::{ReferenceSequence, reference_sequence::index::LinearIndex},
};

use self::{bins::write_bins, intervals::write_intervals, metadata::write_metadata};
use crate::io::writer::num::write_u32_le;

pub(super) fn write_reference_sequences<W>(
    writer: &mut W,
    reference_sequences: &[ReferenceSequence<LinearIndex>],
) -> io::Result<()>
where
    W: Write,
{
    let n_ref = u32::try_from(reference_sequences.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u32_le(writer, n_ref)?;

    for reference_sequence in reference_sequences {
        write_reference_sequence(writer, reference_sequence)?;
    }

    Ok(())
}

fn write_reference_sequence<W>(
    writer: &mut W,
    reference_sequence: &ReferenceSequence<LinearIndex>,
) -> io::Result<()>
where
    W: Write,
{
    write_bins(
        writer,
        reference_sequence.bins(),
        reference_sequence.metadata(),
    )?;

    write_intervals(writer, reference_sequence.index())?;

    Ok(())
}
