mod container;
mod file_id;
mod format_version;
mod magic_number;

use std::io::{self, Write};

use noodles_sam as sam;

use self::{
    container::write_container, file_id::write_file_id, format_version::write_format_version,
    magic_number::write_magic_number,
};
use crate::FileDefinition;

pub fn write_file_header<W>(writer: &mut W, header: &sam::Header) -> io::Result<()>
where
    W: Write,
{
    write_container(writer, header)
}

pub fn write_file_definition<W>(writer: &mut W, file_definition: &FileDefinition) -> io::Result<()>
where
    W: Write,
{
    write_magic_number(writer)?;
    write_format_version(writer, file_definition.version())?;
    write_file_id(writer, file_definition.file_id())?;
    Ok(())
}
