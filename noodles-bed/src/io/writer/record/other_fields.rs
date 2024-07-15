use std::io::{self, Write};

use super::write_separator;
use crate::feature::record::OtherFields;

pub(super) fn write_other_fields<W, F>(writer: &mut W, other_fields: &F) -> io::Result<()>
where
    W: Write,
    F: OtherFields,
{
    for field in other_fields.iter() {
        write_separator(writer)?;
        writer.write_all(field)?;
    }

    Ok(())
}
