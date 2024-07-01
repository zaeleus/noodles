use std::io::{self, Write};

use super::write_separator;
use crate::record::OtherFields;

pub(super) fn write_other_fields<W>(writer: &mut W, other_fields: OtherFields<'_>) -> io::Result<()>
where
    W: Write,
{
    for field in other_fields.iter() {
        write_separator(writer)?;
        writer.write_all(field)?;
    }

    Ok(())
}
