//! BAM record data field value.

pub mod array;

use std::io;

use bstr::{BStr, ByteSlice};
use noodles_sam::alignment::record::data::field::{Type, Value};

use self::array::decode_array;

pub(crate) fn decode_value<'a>(src: &mut &'a [u8], ty: Type) -> io::Result<Value<'a>> {
    match ty {
        Type::Character => read_u8(src).map(Value::Character),
        Type::Int8 => read_u8(src).map(|n| Value::Int8(n as i8)),
        Type::UInt8 => read_u8(src).map(Value::UInt8),
        Type::Int16 => read_u16_le(src).map(|n| Value::Int16(n as i16)),
        Type::UInt16 => read_u16_le(src).map(Value::UInt16),
        Type::Int32 => read_u32_le(src).map(|n| Value::Int32(n as i32)),
        Type::UInt32 => read_u32_le(src).map(Value::UInt32),
        Type::Float => read_f32_le(src).map(Value::Float),
        Type::String => read_string(src).map(Value::String),
        Type::Hex => read_string(src).map(Value::Hex),
        Type::Array => decode_array(src).map(Value::Array),
    }
}

fn read_u8(src: &mut &[u8]) -> io::Result<u8> {
    let Some((n, rest)) = src.split_first() else {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    };

    *src = rest;

    Ok(*n)
}

fn read_u16_le(src: &mut &[u8]) -> io::Result<u16> {
    let Some((buf, rest)) = src.split_first_chunk() else {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    };

    *src = rest;

    Ok(u16::from_le_bytes(*buf))
}

fn read_u32_le(src: &mut &[u8]) -> io::Result<u32> {
    let Some((buf, rest)) = src.split_first_chunk() else {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    };

    *src = rest;

    Ok(u32::from_le_bytes(*buf))
}

fn read_f32_le(src: &mut &[u8]) -> io::Result<f32> {
    let Some((buf, rest)) = src.split_first_chunk() else {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    };

    *src = rest;

    Ok(f32::from_le_bytes(*buf))
}

fn read_string<'a>(src: &mut &'a [u8]) -> io::Result<&'a BStr> {
    const NUL: u8 = 0x00;

    let len = src
        .as_bstr()
        .find_byte(NUL)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "string not NUL terminated"))?;

    let (buf, rest) = src.split_at(len);

    *src = &rest[1..];

    Ok(buf.as_bstr())
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record::data::field::value::Array;

    use super::*;

    #[test]
    fn test_ty() {
        assert_eq!(Value::Character(b'n').ty(), Type::Character);
        assert_eq!(Value::Int8(0).ty(), Type::Int8);
        assert_eq!(Value::UInt8(0).ty(), Type::UInt8);
        assert_eq!(Value::Int16(0).ty(), Type::Int16);
        assert_eq!(Value::UInt16(0).ty(), Type::UInt16);
        assert_eq!(Value::Int32(0).ty(), Type::Int32);
        assert_eq!(Value::UInt32(0).ty(), Type::UInt32);
        assert_eq!(Value::Float(0.0).ty(), Type::Float);
        assert_eq!(Value::String(b"ndls".as_bstr()).ty(), Type::String);
        assert_eq!(Value::Hex(b"CAFE".as_bstr()).ty(), Type::Hex);
        assert_eq!(
            Value::Array(Array::UInt8(Box::new(array::Values::new(&[0])))).ty(),
            Type::Array
        );
    }

    #[test]
    fn test_decode_value() -> io::Result<()> {
        let mut src = &[b'n'][..];
        assert!(matches!(
            decode_value(&mut src, Type::Character)?,
            Value::Character(b'n')
        ));

        let mut src = &[0x00][..];
        assert!(matches!(
            decode_value(&mut src, Type::Int8)?,
            Value::Int8(0)
        ));

        let mut src = &[0x00][..];
        assert!(matches!(
            decode_value(&mut src, Type::UInt8)?,
            Value::UInt8(0)
        ));

        let mut src = &[0x00, 0x00][..];
        assert!(matches!(
            decode_value(&mut src, Type::Int16)?,
            Value::Int16(0)
        ));

        let mut src = &[0x00, 0x00][..];
        assert!(matches!(
            decode_value(&mut src, Type::UInt16)?,
            Value::UInt16(0)
        ));

        let mut src = &[0x00, 0x00, 0x00, 0x00][..];
        assert!(matches!(
            decode_value(&mut src, Type::Int32)?,
            Value::Int32(0)
        ));

        let mut src = &[0x00, 0x00, 0x00, 0x00][..];
        assert!(matches!(
            decode_value(&mut src, Type::UInt32)?,
            Value::UInt32(0)
        ));

        let mut src = &[0x00, 0x00, 0x00, 0x00][..];
        assert!(matches!(
            decode_value(&mut src, Type::Float)?,
            Value::Float(n) if n == 0.0
        ));

        let mut src = &[b'n', b'd', b'l', b's', 0x00][..];
        assert!(matches!(
            decode_value(&mut src, Type::String)?,
            Value::String(s) if s == "ndls"
        ));

        let mut src = &[b'C', b'A', b'F', b'E', 0x00][..];
        assert!(matches!(
            decode_value(&mut src, Type::Hex)?,
            Value::Hex(s) if s == "CAFE"
        ));

        let mut src = &[b'C', 0x01, 0x00, 0x00, 0x00, 0x00][..];
        assert!(matches!(
            decode_value(&mut src, Type::Array)?,
            Value::Array(_)
        ));

        Ok(())
    }
}
