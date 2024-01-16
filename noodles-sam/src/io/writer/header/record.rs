mod kind;
mod value;

use std::io::{self, Write};

use self::kind::write_kind;
use crate::header::record::{
    value::{
        map::{Header, Program, ReadGroup, ReferenceSequence},
        Map,
    },
    Kind,
};

const DELIMITER: u8 = b'\t';
const LINE_FEED: u8 = b'\n';
const PREFIX: u8 = b'@';

pub(super) fn write_header<W>(writer: &mut W, header: &Map<Header>) -> io::Result<()>
where
    W: Write,
{
    write_prefix(writer)?;
    write_kind(writer, Kind::Header)?;
    value::map::write_header(writer, header)?;
    write_newline(writer)?;
    Ok(())
}

pub(super) fn write_reference_sequence<W>(
    writer: &mut W,
    name: &[u8],
    reference_sequence: &Map<ReferenceSequence>,
) -> io::Result<()>
where
    W: Write,
{
    write_prefix(writer)?;
    write_kind(writer, Kind::ReferenceSequence)?;
    value::map::write_reference_sequence(writer, name, reference_sequence)?;
    write_newline(writer)?;
    Ok(())
}

pub(super) fn write_read_group<W>(
    writer: &mut W,
    id: &[u8],
    read_group: &Map<ReadGroup>,
) -> io::Result<()>
where
    W: Write,
{
    write_prefix(writer)?;
    write_kind(writer, Kind::ReadGroup)?;
    value::map::write_read_group(writer, id, read_group)?;
    write_newline(writer)?;
    Ok(())
}

pub(super) fn write_program<W>(writer: &mut W, id: &str, program: &Map<Program>) -> io::Result<()>
where
    W: Write,
{
    write_prefix(writer)?;
    write_kind(writer, Kind::Program)?;
    value::map::write_program(writer, id, program)?;
    write_newline(writer)?;
    Ok(())
}

pub(super) fn write_comment<W>(writer: &mut W, comment: &str) -> io::Result<()>
where
    W: Write,
{
    write_prefix(writer)?;
    write_kind(writer, Kind::Comment)?;
    write_delimiter(writer)?;
    value::write_string(writer, comment)?;
    write_newline(writer)?;
    Ok(())
}

fn write_prefix<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&[PREFIX])
}

fn write_delimiter<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&[DELIMITER])
}

fn write_newline<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&[LINE_FEED])
}
