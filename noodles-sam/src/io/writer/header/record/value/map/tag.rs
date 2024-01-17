use std::io::{self, Write};

pub(super) fn write_tag<W>(writer: &mut W, tag: [u8; 2]) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&tag)
}
