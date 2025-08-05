mod itf8;
mod ltf8;
mod vlq;

use std::io::{self, Write};

pub use self::{itf8::write_itf8, ltf8::write_ltf8, vlq::write_uint7};

pub(crate) fn write_u8<W>(writer: &mut W, n: u8) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&[n])
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
