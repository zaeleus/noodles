use std::io::{self, Write};

use crate::fai::Record;

pub(crate) fn write_record<W>(writer: &mut W, record: &Record) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(record.name())?;

    writeln!(
        writer,
        "\t{length}\t{offset}\t{line_base_count}\t{line_width}",
        length = record.length(),
        offset = record.position(),
        line_base_count = record.line_base_count(),
        line_width = record.line_width(),
    )
}
