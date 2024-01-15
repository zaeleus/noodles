mod sequence;

use std::io::{self, Write};

use self::sequence::write_sequence;
use crate::Record;

pub(super) fn write_record<W>(
    writer: &mut W,
    record: &Record,
    line_base_count: usize,
) -> io::Result<()>
where
    W: Write,
{
    writeln!(writer, "{}", record.definition())?;
    write_sequence(writer, record.sequence(), line_base_count)?;
    Ok(())
}
