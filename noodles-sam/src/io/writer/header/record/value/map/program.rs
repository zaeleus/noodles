use std::io::{self, Write};

use super::{fmt_display_field, write_other_fields};
use crate::header::record::value::{
    map::{program::tag, Program},
    Map,
};

pub(crate) fn write_program<W>(writer: &mut W, id: &[u8], program: &Map<Program>) -> io::Result<()>
where
    W: Write,
{
    use crate::io::writer::header::record::{value::map::write_separator, write_delimiter};

    write_delimiter(writer)?;
    writer.write_all(tag::ID.as_ref())?;
    write_separator(writer)?;
    writer.write_all(id)?;

    if let Some(name) = program.name() {
        fmt_display_field(writer, tag::NAME, name)?;
    }

    if let Some(command_line) = program.command_line() {
        fmt_display_field(writer, tag::COMMAND_LINE, command_line)?;
    }

    if let Some(previous_id) = program.previous_id() {
        fmt_display_field(writer, tag::PREVIOUS_ID, previous_id)?;
    }

    if let Some(description) = program.description() {
        fmt_display_field(writer, tag::DESCRIPTION, description)?;
    }

    if let Some(version) = program.version() {
        fmt_display_field(writer, tag::VERSION, version)?;
    }

    write_other_fields(writer, program.other_fields())?;

    Ok(())
}
