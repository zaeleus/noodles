mod length;
mod name;

use std::io::{self, Write};

use self::{length::write_length_field, name::write_name_field};
use super::write_other_fields;
use crate::header::record::value::{map::ReferenceSequence, Map};

pub(crate) fn write_reference_sequence<W>(
    writer: &mut W,
    name: &[u8],
    reference_sequence: &Map<ReferenceSequence>,
) -> io::Result<()>
where
    W: Write,
{
    write_name_field(writer, name)?;
    write_length_field(writer, reference_sequence.length())?;
    write_other_fields(writer, reference_sequence.other_fields())?;
    Ok(())
}
