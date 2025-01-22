mod bins;
mod intervals;
mod metadata;

use std::io::{self, Write};

use noodles_csi::binning_index::{
    index::{reference_sequence::index::LinearIndex, ReferenceSequence},
    ReferenceSequence as _,
};

use self::{bins::write_bins, intervals::write_intervals, metadata::write_metadata};

pub(super) fn write_reference_sequences<W>(
    writer: &mut W,
    reference_sequences: &[ReferenceSequence<LinearIndex>],
) -> io::Result<()>
where
    W: Write,
{
    for reference_sequence in reference_sequences {
        write_reference_sequence(writer, reference_sequence)?;
    }

    Ok(())
}

pub fn write_reference_sequence<W>(
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
