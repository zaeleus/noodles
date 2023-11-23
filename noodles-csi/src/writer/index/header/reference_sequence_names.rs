use std::io::{self, Write};

use byteorder::{LittleEndian, WriteBytesExt};

use crate::binning_index::index::header::ReferenceSequenceNames;

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
    writer.write_i32::<LittleEndian>(l_nm)?;

    for reference_sequence_name in reference_sequence_names {
        writer.write_all(reference_sequence_name.as_bytes())?;
        writer.write_u8(NUL)?;
    }

    Ok(())
}
