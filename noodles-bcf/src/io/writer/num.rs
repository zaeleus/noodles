use std::io::{self, Write};

pub(crate) fn write_i8<W>(writer: &mut W, n: i8) -> io::Result<()>
where
    W: Write,
{
    write_u8(writer, n as u8)
}

pub(crate) fn write_u8<W>(writer: &mut W, n: u8) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&[n])
}

pub(crate) fn write_i16_le<W>(writer: &mut W, n: i16) -> io::Result<()>
where
    W: Write,
{
    let buf = n.to_le_bytes();
    writer.write_all(&buf)
}

pub(crate) fn write_u16_le<W>(writer: &mut W, n: u16) -> io::Result<()>
where
    W: Write,
{
    let buf = n.to_le_bytes();
    writer.write_all(&buf)
}

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

pub(crate) fn write_f32_le<W>(writer: &mut W, n: f32) -> io::Result<()>
where
    W: Write,
{
    let buf = n.to_le_bytes();
    writer.write_all(&buf)
}
