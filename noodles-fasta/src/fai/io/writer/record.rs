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

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use super::*;

    #[test]
    fn test_write_record() -> io::Result<()> {
        let mut buf = Vec::new();

        let line_base_count = const { NonZero::new(80).unwrap() };
        let line_width = const { NonZero::new(81).unwrap() };
        let record = Record::new("sq0", 13, 5, line_base_count, line_width);

        write_record(&mut buf, &record)?;

        assert_eq!(buf, b"sq0\t13\t5\t80\t81\n");

        Ok(())
    }
}
