mod itf8;
mod ltf8;
mod vlq;

use std::io::{self, Write};

pub use self::{itf8::write_itf8, ltf8::write_ltf8, vlq::write_uint7};

pub(crate) fn write_i32_le<W>(writer: &mut W, n: i32) -> io::Result<()>
where
    W: Write,
{
    let buf = n.to_le_bytes();
    writer.write_all(&buf)
}
