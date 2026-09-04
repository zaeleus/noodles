use std::{
    io::{self, Read, Write},
    mem,
};

pub(crate) fn read_u32_le<R>(reader: &mut R) -> io::Result<u32>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u32>()];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

pub(crate) fn read_u64_le<R>(reader: &mut R) -> io::Result<u64>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u64>()];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
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
