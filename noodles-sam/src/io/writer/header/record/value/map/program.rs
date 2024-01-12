use std::io::{self, Write};

use super::{write_field, write_other_fields};
use crate::header::record::value::{
    map::{program::tag, Program},
    Map,
};

pub(crate) fn write_program<W>(writer: &mut W, id: &str, program: &Map<Program>) -> io::Result<()>
where
    W: Write,
{
    write_field(writer, tag::ID, id)?;

    if let Some(name) = program.name() {
        write_field(writer, tag::NAME, name)?;
    }

    if let Some(command_line) = program.command_line() {
        write_field(writer, tag::COMMAND_LINE, command_line)?;
    }

    if let Some(previous_id) = program.previous_id() {
        write_field(writer, tag::PREVIOUS_ID, previous_id)?;
    }

    if let Some(description) = program.description() {
        write_field(writer, tag::DESCRIPTION, description)?;
    }

    if let Some(version) = program.version() {
        write_field(writer, tag::VERSION, version)?;
    }

    write_other_fields(writer, program.other_fields())?;

    Ok(())
}
