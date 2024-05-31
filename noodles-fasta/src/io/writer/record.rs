mod definition;
mod sequence;

use std::io::{self, Write};

use self::{definition::write_definition, sequence::write_sequence};
use crate::Record;

pub(super) fn write_record<W>(
    writer: &mut W,
    record: &Record,
    line_base_count: usize,
) -> io::Result<()>
where
    W: Write,
{
    write_definition(writer, record.definition())?;
    write_newline(writer)?;

    write_sequence(writer, record.sequence(), line_base_count)?;

    Ok(())
}

fn write_newline<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const LINE_FEED: u8 = b'\n';
    writer.write_all(&[LINE_FEED])
}
