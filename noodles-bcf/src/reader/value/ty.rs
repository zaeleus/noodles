use std::io::{self, Read};

use byteorder::ReadBytesExt;

use crate::lazy::record::{value::Type, Value};

use super::read_value;

pub fn read_type<R>(reader: &mut R) -> io::Result<Option<Type>>
where
    R: Read,
{
    use super::{Int16, Int32, Int8};

    let encoding = reader.read_u8()?;

    let mut len = usize::from(encoding >> 4);

    if len == 0x0f {
        let value = read_value(reader)?;

        let next_len = match value {
            Some(Value::Int8(Some(Int8::Value(n)))) => i32::from(n),
            Some(Value::Int16(Some(Int16::Value(n)))) => i32::from(n),
            Some(Value::Int32(Some(Int32::Value(n)))) => n,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid length value: #{value:?}"),
                ))
            }
        };

        len =
            usize::try_from(next_len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    let ty = encoding & 0x0f;

    match ty {
        0 => Ok(None),
        1 => Ok(Some(Type::Int8(len))),
        2 => Ok(Some(Type::Int16(len))),
        3 => Ok(Some(Type::Int32(len))),
        5 => Ok(Some(Type::Float(len))),
        7 => Ok(Some(Type::String(len))),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid type: {ty}"),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_type() -> io::Result<()> {
        let data = [0x00];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, None);

        let data = [0x11];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Some(Type::Int8(1)));

        let data = [0x31];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Some(Type::Int8(3)));

        let data = [0x12];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Some(Type::Int16(1)));

        let data = [0x32];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Some(Type::Int16(3)));

        let data = [0x13];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Some(Type::Int32(1)));

        let data = [0x33];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Some(Type::Int32(3)));

        let data = [0x15];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Some(Type::Float(1)));

        let data = [0x35];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Some(Type::Float(3)));

        let data = [0x17];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Some(Type::String(1)));

        let data = [0x37];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Some(Type::String(3)));

        Ok(())
    }

    #[test]
    fn test_read_type_with_gt_14_values() -> io::Result<()> {
        let data = [0xf1, 0x11, 0x15];
        let mut reader = &data[..];
        assert_eq!(read_type(&mut reader)?, Some(Type::Int8(21)));
        Ok(())
    }
}
