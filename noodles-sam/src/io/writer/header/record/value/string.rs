use std::io::{self, Write};

pub(crate) fn write_string<W>(writer: &mut W, buf: &[u8]) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(buf)
}
