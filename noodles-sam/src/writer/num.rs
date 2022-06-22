use std::io::{self, Write};

pub(super) fn write_int<W, I>(writer: &mut W, i: I) -> io::Result<()>
where
    W: Write,
    I: itoa::Integer,
{
    let mut buf = itoa::Buffer::new();
    let a = buf.format(i);
    writer.write_all(a.as_bytes())
}
