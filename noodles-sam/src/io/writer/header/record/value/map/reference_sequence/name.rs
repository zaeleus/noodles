use std::io::{self, Write};

use crate::header::record::value::map::reference_sequence::tag;

pub(super) fn write_name_field<W>(writer: &mut W, name: &[u8]) -> io::Result<()>
where
    W: Write,
{
    use crate::io::writer::header::record::{value::map::write_separator, write_delimiter};

    write_delimiter(writer)?;
    writer.write_all(tag::NAME.as_ref())?;
    write_separator(writer)?;
    writer.write_all(name)?;

    Ok(())
}
