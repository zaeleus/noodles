use std::io::{self, Write};

use crate::Record;

pub(super) fn write_record<W>(
    writer: &mut W,
    definition_separator: u8,
    record: &Record,
) -> io::Result<()>
where
    W: Write,
{
    const NAME_PREFIX: &[u8] = b"@";
    const LINE_FEED: &[u8] = b"\n";

    writer.write_all(NAME_PREFIX)?;
    writer.write_all(record.name())?;

    if !record.description().is_empty() {
        writer.write_all(&[definition_separator])?;
        writer.write_all(record.description())?;
    }

    writer.write_all(LINE_FEED)?;

    writer.write_all(record.sequence())?;
    writer.write_all(LINE_FEED)?;

    writer.write_all(b"+")?;
    writer.write_all(LINE_FEED)?;

    writer.write_all(record.quality_scores())?;
    writer.write_all(LINE_FEED)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record() -> io::Result<()> {
        use crate::record::Definition;

        const HORIZONTAL_TAB: u8 = b'\t';
        const SPACE: u8 = b' ';

        let mut record = Record::new(Definition::new("r0", ""), "ACGT", "NDLS");

        let mut buf = Vec::new();
        write_record(&mut buf, SPACE, &record)?;
        let expected = b"@r0\nACGT\n+\nNDLS\n";
        assert_eq!(buf, expected);

        record.description_mut().extend_from_slice(b"LN:4");

        buf.clear();
        write_record(&mut buf, SPACE, &record)?;
        let expected = b"@r0 LN:4\nACGT\n+\nNDLS\n";
        assert_eq!(buf, expected);

        buf.clear();
        write_record(&mut buf, HORIZONTAL_TAB, &record)?;
        let expected = b"@r0\tLN:4\nACGT\n+\nNDLS\n";
        assert_eq!(buf, expected);

        Ok(())
    }
}
