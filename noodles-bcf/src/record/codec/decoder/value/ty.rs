use std::io;

use byteorder::ReadBytesExt;

use crate::lazy::record::value::Type;

use super::read_value;

pub fn read_type(src: &mut &[u8]) -> io::Result<Option<Type>> {
    let encoding = src.read_u8()?;

    let mut len = usize::from(encoding >> 4);

    if len == 0x0f {
        let value = read_value(src)?;

        len = match value.as_ref().and_then(|v| v.as_int()) {
            Some(n) => {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid length value: #{value:?}"),
                ))
            }
        };
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
        let mut src = &[0x00][..];
        assert_eq!(read_type(&mut src)?, None);

        let mut src = &[0x11][..];
        assert_eq!(read_type(&mut src)?, Some(Type::Int8(1)));

        let mut src = &[0x31][..];
        assert_eq!(read_type(&mut src)?, Some(Type::Int8(3)));

        let mut src = &[0x12][..];
        assert_eq!(read_type(&mut src)?, Some(Type::Int16(1)));

        let mut src = &[0x32][..];
        assert_eq!(read_type(&mut src)?, Some(Type::Int16(3)));

        let mut src = &[0x13][..];
        assert_eq!(read_type(&mut src)?, Some(Type::Int32(1)));

        let mut src = &[0x33][..];
        assert_eq!(read_type(&mut src)?, Some(Type::Int32(3)));

        let mut src = &[0x15][..];
        assert_eq!(read_type(&mut src)?, Some(Type::Float(1)));

        let mut src = &[0x35][..];
        assert_eq!(read_type(&mut src)?, Some(Type::Float(3)));

        let mut src = &[0x17][..];
        assert_eq!(read_type(&mut src)?, Some(Type::String(1)));

        let mut src = &[0x37][..];
        assert_eq!(read_type(&mut src)?, Some(Type::String(3)));

        Ok(())
    }

    #[test]
    fn test_read_type_with_gt_14_values() -> io::Result<()> {
        let mut src = &[0xf1, 0x11, 0x15][..];
        assert_eq!(read_type(&mut src)?, Some(Type::Int8(21)));
        Ok(())
    }
}
