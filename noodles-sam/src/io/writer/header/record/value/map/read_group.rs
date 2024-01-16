use std::io::{self, Write};

use super::{write_field, write_other_fields};
use crate::header::record::value::{
    map::{read_group::tag, ReadGroup},
    Map,
};

pub(crate) fn write_read_group<W>(
    writer: &mut W,
    id: &[u8],
    read_group: &Map<ReadGroup>,
) -> io::Result<()>
where
    W: Write,
{
    write_field(writer, tag::ID, id)?;
    write_other_fields(writer, read_group.other_fields())?;
    Ok(())
}
