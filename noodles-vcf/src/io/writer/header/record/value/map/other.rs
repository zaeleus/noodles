use std::io::{self, Write};

use super::write_other_fields;
use crate::header::record::value::{map::Other, Map};

pub(crate) fn write_other<W>(writer: &mut W, other: &Map<Other>) -> io::Result<()>
where
    W: Write,
{
    write_other_fields(writer, other.other_fields())
}
