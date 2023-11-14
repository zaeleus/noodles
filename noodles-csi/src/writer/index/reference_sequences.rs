mod bins;
mod metadata;

use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use self::{bins::write_bins, metadata::write_metadata};
use crate::index::ReferenceSequence;

pub(super) fn write_reference_sequences<W>(
    writer: &mut W,
    depth: u8,
    reference_sequences: &[ReferenceSequence],
) -> io::Result<()>
where
    W: Write,
{
    let n_ref = i32::try_from(reference_sequences.len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    writer.write_i32::<LittleEndian>(n_ref)?;

    for reference_sequence in reference_sequences {
        write_bins(
            writer,
            depth,
            reference_sequence.bins(),
            reference_sequence.metadata(),
        )?;
    }

    Ok(())
}
