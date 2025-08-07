mod bins;
mod metadata;

use std::io::{self, Write};

use self::{bins::write_bins, metadata::write_metadata};
use crate::{
    binning_index::{
        ReferenceSequence as _,
        index::{ReferenceSequence, reference_sequence::index::BinnedIndex},
    },
    io::writer::num::write_i32_le,
};

pub(super) fn write_reference_sequences<W>(
    writer: &mut W,
    depth: u8,
    reference_sequences: &[ReferenceSequence<BinnedIndex>],
) -> io::Result<()>
where
    W: Write,
{
    let n_ref = i32::try_from(reference_sequences.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, n_ref)?;

    for reference_sequence in reference_sequences {
        write_bins(
            writer,
            depth,
            reference_sequence.bins(),
            reference_sequence.index(),
            reference_sequence.metadata(),
        )?;
    }

    Ok(())
}
