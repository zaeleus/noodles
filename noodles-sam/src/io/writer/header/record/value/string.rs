use std::io::{self, Write};

pub(crate) fn write_string<W>(writer: &mut W, s: &str) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(s.as_bytes())
}
