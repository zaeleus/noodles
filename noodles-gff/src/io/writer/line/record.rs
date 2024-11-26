mod phase;
mod position;
mod score;
mod strand;

use std::io::{self, Write};

use self::{
    phase::write_phase, position::write_position, score::write_score, strand::write_strand,
};
use crate::RecordBuf;

pub(crate) fn write_record<W>(writer: &mut W, record: &RecordBuf) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(record.reference_sequence_name().as_bytes())?;
    write_separator(writer)?;

    writer.write_all(record.source().as_bytes())?;
    write_separator(writer)?;

    writer.write_all(record.ty().as_bytes())?;
    write_separator(writer)?;

    write_position(writer, record.start())?;
    write_separator(writer)?;

    write_position(writer, record.end())?;
    write_separator(writer)?;

    write_score(writer, record.score())?;
    write_separator(writer)?;

    write_strand(writer, record.strand())?;
    write_separator(writer)?;

    write_phase(writer, record.phase())?;
    write_separator(writer)?;

    let attributes = record.attributes();

    if attributes.is_empty() {
        write_missing(writer)?;
    } else {
        write!(writer, "{}", attributes)?;
    }

    Ok(())
}

fn write_missing<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const MISSING: u8 = b'.';
    writer.write_all(&[MISSING])
}

fn write_separator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const SEPARATOR: u8 = b'\t';
    writer.write_all(&[SEPARATOR])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record() -> io::Result<()> {
        let mut buf = Vec::new();
        let record = RecordBuf::default();
        write_record(&mut buf, &record)?;
        assert_eq!(buf, b".\t.\t.\t1\t1\t.\t.\t.\t.");
        Ok(())
    }
}
