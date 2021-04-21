use std::io::{self, Read};

use byteorder::ReadBytesExt;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Type {
    Int8,
    Int8Array(usize),
    Int16,
    Int16Array(usize),
    Int32,
    Int32Array(usize),
    Float,
    FloatArray(usize),
    Char,
    String(usize),
}

pub fn read_type<R>(reader: &mut R) -> io::Result<Type>
where
    R: Read,
{
    let encoding = reader.read_u8()?;

    let len = usize::from(encoding >> 4);
    assert!(len < 0x0f);

    let ty = encoding & 0x0f;

    match ty {
        1 => {
            if len == 1 {
                Ok(Type::Int8)
            } else {
                Ok(Type::Int8Array(len))
            }
        }
        2 => {
            if len == 1 {
                Ok(Type::Int16)
            } else {
                Ok(Type::Int16Array(len))
            }
        }
        3 => {
            if len == 1 {
                Ok(Type::Int32)
            } else {
                Ok(Type::Int32Array(len))
            }
        }
        5 => {
            if len == 1 {
                Ok(Type::Float)
            } else {
                Ok(Type::FloatArray(len))
            }
        }
        7 => {
            if len == 1 {
                Ok(Type::Char)
            } else {
                Ok(Type::String(len))
            }
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid type: {}", ty),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_type() -> io::Result<()> {
        let data = [0x11];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Int8);

        let data = [0x31];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Int8Array(3));

        let data = [0x12];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Int16);

        let data = [0x32];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Int16Array(3));

        let data = [0x13];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Int32);

        let data = [0x33];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Int32Array(3));

        let data = [0x15];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Float);

        let data = [0x35];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::FloatArray(3));

        let data = [0x17];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Char);

        let data = [0x37];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::String(3));

        Ok(())
    }
}
