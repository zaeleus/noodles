use std::{
    convert::TryFrom,
    io::{self, Read},
};

use byteorder::ReadBytesExt;

use super::{read_value, Value};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Type {
    Int8(usize),
    Int16(usize),
    Int32(usize),
    Float(usize),
    String(usize),
}

pub fn read_type<R>(reader: &mut R) -> io::Result<Type>
where
    R: Read,
{
    let encoding = reader.read_u8()?;

    let mut len = usize::from(encoding >> 4);

    if len == 0x0f {
        let value = read_value(reader)?;

        let next_len = match value {
            Value::Int8(Some(n)) => i32::from(n),
            Value::Int16(Some(n)) => i32::from(n),
            Value::Int32(Some(n)) => n,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid length value: #{:?}", value),
                ))
            }
        };

        len =
            usize::try_from(next_len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    let ty = encoding & 0x0f;

    match ty {
        1 => Ok(Type::Int8(len)),
        2 => Ok(Type::Int16(len)),
        3 => Ok(Type::Int32(len)),
        5 => Ok(Type::Float(len)),
        7 => Ok(Type::String(len)),
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
        assert_eq!(read_type(&mut reader)?, Type::Int8(1));

        let data = [0x31];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Int8(3));

        let data = [0x12];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Int16(1));

        let data = [0x32];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Int16(3));

        let data = [0x13];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Int32(1));

        let data = [0x33];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Int32(3));

        let data = [0x15];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Float(1));

        let data = [0x35];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Float(3));

        let data = [0x17];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::String(1));

        let data = [0x37];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::String(3));

        Ok(())
    }

    #[test]
    fn test_read_type_with_gt_14_values() -> io::Result<()> {
        let data = [0xf1, 0x11, 0x15];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Type::Int8(21));
        Ok(())
    }
}
