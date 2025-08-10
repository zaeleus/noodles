use std::io::{self, Write};

use byteorder::WriteBytesExt;
use noodles_csi::binning_index::index::header::ReferenceSequenceNames;

use crate::io::writer::num::write_i32_le;

pub(super) fn write_reference_sequence_names<W>(
    writer: &mut W,
    reference_sequence_names: &ReferenceSequenceNames,
) -> io::Result<()>
where
    W: Write,
{
    const NUL: u8 = 0x00;

    // Add 1 for each trailing NUL.
    let len = reference_sequence_names
        .iter()
        .map(|n| n.len() + 1)
        .sum::<usize>();
    let l_nm = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_i32_le(writer, l_nm)?;

    for reference_sequence_name in reference_sequence_names {
        writer.write_all(reference_sequence_name)?;
        writer.write_u8(NUL)?;
    }

    Ok(())
}
