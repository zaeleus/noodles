mod header;
mod program;
mod read_group;
mod reference_sequence;
mod tag;

use std::io::{self, Write};

use self::tag::write_tag;
pub(crate) use self::{
    header::write_header, program::write_program, read_group::write_read_group,
    reference_sequence::write_reference_sequence,
};
use crate::header::record::value::map::OtherFields;

const SEPARATOR: u8 = b':';

fn write_field<W, K>(writer: &mut W, tag: K, value: &[u8]) -> io::Result<()>
where
    W: Write,
    K: AsRef<[u8; 2]>,
{
    use crate::io::writer::header::record::write_delimiter;

    write_delimiter(writer)?;
    write_tag(writer, *tag.as_ref())?;
    write_separator(writer)?;
    writer.write_all(value)?;

    Ok(())
}

fn write_other_fields<W, S>(writer: &mut W, other_fields: &OtherFields<S>) -> io::Result<()>
where
    W: Write,
{
    for (tag, value) in other_fields {
        write_field(writer, tag, value)?;
    }

    Ok(())
}

fn write_separator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&[SEPARATOR])
}
