mod value;

use std::io::{self, Write};

use self::value::write_value;
use super::write_separator;
use crate::feature::record::OtherFields;

pub(super) fn write_other_fields<W, F>(writer: &mut W, other_fields: &F) -> io::Result<()>
where
    W: Write,
    F: OtherFields + ?Sized,
{
    for value in other_fields.iter() {
        write_separator(writer)?;
        write_value(writer, value)?;
    }

    Ok(())
}
