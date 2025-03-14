mod attributes;
mod frame;
mod position;
mod score;
mod strand;

use std::io::{self, Write};

use self::{
    attributes::write_attributes, frame::write_frame, position::write_position, score::write_score,
    strand::write_strand,
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
    write_frame(writer, record.frame())?;

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
