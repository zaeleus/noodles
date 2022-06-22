use std::io::{self, Write};

use lexical_core::FormattedSize;

pub fn write_i8<W>(writer: &mut W, n: i8) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; i8::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_u8<W>(writer: &mut W, n: u8) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; u8::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_i16<W>(writer: &mut W, n: i16) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; i16::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_u16<W>(writer: &mut W, n: u16) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; u16::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_i32<W>(writer: &mut W, n: i32) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; i32::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_u32<W>(writer: &mut W, n: u32) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; u32::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}

pub fn write_usize<W>(writer: &mut W, n: usize) -> io::Result<()>
where
    W: Write,
{
    let mut dst = [0; usize::FORMATTED_SIZE_DECIMAL];
    let buf = lexical_core::write(n, &mut dst);
    writer.write_all(buf)
}
