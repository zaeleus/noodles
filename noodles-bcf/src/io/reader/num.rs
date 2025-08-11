use std::{
    io::{self, Read},
    mem,
};

pub(crate) fn read_u16_le<R>(reader: &mut R) -> io::Result<u16>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u16>()];
    reader.read_exact(&mut buf)?;
    Ok(u16::from_le_bytes(buf))
}

pub(crate) fn read_i32_le<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<i32>()];
    reader.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

pub(crate) fn read_u32_le<R>(reader: &mut R) -> io::Result<u32>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u32>()];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

pub(crate) fn read_f32_le<R>(reader: &mut R) -> io::Result<f32>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<f32>()];
    reader.read_exact(&mut buf)?;
    Ok(f32::from_le_bytes(buf))
}
