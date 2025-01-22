use std::io::{self, Write};

use crate::MAGIC_NUMBER;

pub(super) fn write_magic_number<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&MAGIC_NUMBER)
}
