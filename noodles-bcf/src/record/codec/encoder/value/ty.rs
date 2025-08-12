use std::io::{self, Write};

use super::write_value;
use crate::{
    io::writer::num::write_u8,
    record::codec::value::{Int8, Int16, Int32, Type, Value},
};

// § 6.3.3 Type encoding (2021-01-13)
const MAX_TYPE_LEN: usize = 0x0f;

pub fn write_type<W>(writer: &mut W, ty: Option<Type>) -> io::Result<()>
where
    W: Write,
{
    let (raw_ty, ty_len): (u8, usize) = match ty {
        None => (0, 0),
        Some(Type::Int8(len)) => (1, len),
        Some(Type::Int16(len)) => (2, len),
        Some(Type::Int32(len)) => (3, len),
        Some(Type::Float(len)) => (5, len),
        Some(Type::String(len)) => (7, len),
    };

    let len = ty_len.clamp(0, MAX_TYPE_LEN) as u8;
    let encoding = (len << 4) | raw_ty;

    write_u8(writer, encoding)?;

    if ty_len >= MAX_TYPE_LEN {
        let len_value = if let Ok(n) = i8::try_from(ty_len) {
            Value::Int8(Some(Int8::Value(n)))
        } else if let Ok(n) = i16::try_from(ty_len) {
            Value::Int16(Some(Int16::Value(n)))
        } else {
            i32::try_from(ty_len)
                .map(|n| Value::Int32(Some(Int32::Value(n))))
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
        };

        write_value(writer, Some(len_value))?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_type() -> io::Result<()> {
        let mut buf = Vec::new();
        write_type(&mut buf, None)?;
        assert_eq!(buf, [0x00]);

        buf.clear();
        write_type(&mut buf, Some(Type::Int8(1)))?;
        assert_eq!(buf, [0x11]);

        buf.clear();
        write_type(&mut buf, Some(Type::Int8(3)))?;
        assert_eq!(buf, [0x31]);

        buf.clear();
        write_type(&mut buf, Some(Type::Int16(1)))?;
        assert_eq!(buf, [0x12]);

        buf.clear();
        write_type(&mut buf, Some(Type::Int16(3)))?;
        assert_eq!(buf, [0x32]);

        buf.clear();
        write_type(&mut buf, Some(Type::Int32(1)))?;
        assert_eq!(buf, [0x13]);

        buf.clear();
        write_type(&mut buf, Some(Type::Int32(3)))?;
        assert_eq!(buf, [0x33]);

        buf.clear();
        write_type(&mut buf, Some(Type::Float(1)))?;
        assert_eq!(buf, [0x15]);

        buf.clear();
        write_type(&mut buf, Some(Type::Float(3)))?;
        assert_eq!(buf, [0x35]);

        buf.clear();
        write_type(&mut buf, Some(Type::String(1)))?;
        assert_eq!(buf, [0x17]);

        buf.clear();
        write_type(&mut buf, Some(Type::String(3)))?;
        assert_eq!(buf, [0x37]);

        Ok(())
    }

    #[test]
    fn test_write_type_with_gt_14_values() -> io::Result<()> {
        let mut buf = Vec::new();
        write_type(&mut buf, Some(Type::Int8(21)))?;
        assert_eq!(buf, [0xf1, 0x11, 0x15]);
        Ok(())
    }
}
