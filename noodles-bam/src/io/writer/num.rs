use std::io::{self, Write};

pub(crate) fn write_i32_le<W>(writer: &mut W, n: i32) -> io::Result<()>
where
    W: Write,
{
    let buf = n.to_le_bytes();
    writer.write_all(&buf)
}

pub(crate) fn write_u32_le<W>(writer: &mut W, n: u32) -> io::Result<()>
where
    W: Write,
{
    let buf = n.to_le_bytes();
    writer.write_all(&buf)
}

pub(crate) fn write_u64_le<W>(writer: &mut W, n: u64) -> io::Result<()>
where
    W: Write,
{
    let buf = n.to_le_bytes();
    writer.write_all(&buf)
}
