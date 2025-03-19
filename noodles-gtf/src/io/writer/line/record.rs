mod attributes;
mod phase;
mod position;
mod reference_sequence_name;
mod score;
mod strand;

use std::io::{self, Write};

use noodles_gff::feature::Record;

use self::{
    attributes::write_attributes, phase::write_phase, position::write_position,
    reference_sequence_name::write_reference_sequence_name, score::write_score,
    strand::write_strand,
};

pub(crate) fn write_record<W, R>(writer: &mut W, record: &R) -> io::Result<()>
where
    W: Write,
    R: Record,
{
    write_reference_sequence_name(writer, record.reference_sequence_name())?;

    write_separator(writer)?;
    writer.write_all(record.source().as_bytes())?;

    write_separator(writer)?;
    writer.write_all(record.ty().as_bytes())?;

    write_separator(writer)?;
    write_position(writer, record.feature_start()?)?;

    write_separator(writer)?;
    write_position(writer, record.feature_end()?)?;

    write_separator(writer)?;
    write_score(writer, record.score().transpose()?)?;

    write_separator(writer)?;
    write_strand(writer, record.strand()?)?;

    write_separator(writer)?;
    write_phase(writer, record.phase().transpose()?)?;

    write_separator(writer)?;
    write_attributes(writer, record.attributes().as_ref())?;

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
