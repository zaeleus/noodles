use std::io::{self, Write};

use super::{write_field, write_other_fields};
use crate::header::record::value::{
    map::{program::tag, Program},
    Map,
};

pub(crate) fn write_program<W>(writer: &mut W, id: &[u8], program: &Map<Program>) -> io::Result<()>
where
    W: Write,
{
    write_field(writer, tag::ID, id)?;
    write_other_fields(writer, program.other_fields())?;
    Ok(())
}
