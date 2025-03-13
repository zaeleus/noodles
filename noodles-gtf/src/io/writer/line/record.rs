mod attributes;
mod position;
mod score;

use std::io::{self, Write};

use self::{attributes::write_attributes, position::write_position, score::write_score};
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

    if let Some(strand) = record.strand() {
        write!(writer, "{strand}")?;
    } else {
        write_missing(writer)?;
    }

    write_separator(writer)?;

    if let Some(frame) = record.frame() {
        write!(writer, "{frame}")?;
    } else {
        write_missing(writer)?;
    }

    write_separator(writer)?;
    write_attributes(writer, record.attributes())?;

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
