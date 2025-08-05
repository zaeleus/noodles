mod itf8;
mod ltf8;
mod vlq;

use std::{
    io::{self, Read},
    mem,
};

pub use self::{
    itf8::{read_itf8, read_itf8_as},
    ltf8::{read_ltf8, read_ltf8_as},
    vlq::read_uint7,
};

pub(crate) fn read_u8<R>(reader: &mut R) -> io::Result<u8>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u8>()];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn read_u16_be<R>(reader: &mut R) -> io::Result<u16>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u16>()];
    reader.read_exact(&mut buf)?;
    Ok(u16::from_be_bytes(buf))
}

pub(crate) fn read_u16_le<R>(reader: &mut R) -> io::Result<u16>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u16>()];
    reader.read_exact(&mut buf)?;
    Ok(u16::from_le_bytes(buf))
}

fn read_u24_be<R>(reader: &mut R) -> io::Result<u32>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u32>()];
    reader.read_exact(&mut buf[1..])?;
    Ok(u32::from_be_bytes(buf))
}

pub(crate) fn read_i32_le<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<i32>()];
    reader.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

fn read_u32_be<R>(reader: &mut R) -> io::Result<u32>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u32>()];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_be_bytes(buf))
}

pub(crate) fn read_u32_le<R>(reader: &mut R) -> io::Result<u32>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u32>()];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}
