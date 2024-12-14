//! BAM record data field value.

pub mod array;

use std::io;

use bstr::{BStr, ByteSlice};
use byteorder::{LittleEndian, ReadBytesExt};
use noodles_sam::alignment::record::data::field::{Type, Value};

use self::array::decode_array;

pub(crate) fn decode_value<'a>(src: &mut &'a [u8], ty: Type) -> io::Result<Value<'a>> {
    match ty {
        Type::Character => read_u8(src).map(Value::Character),
        Type::Int8 => read_u8(src).map(|n| Value::Int8(n as i8)),
        Type::UInt8 => read_u8(src).map(Value::UInt8),
        Type::Int16 => decode_i16(src),
        Type::UInt16 => decode_u16(src),
        Type::Int32 => decode_i32(src),
        Type::UInt32 => decode_u32(src),
        Type::Float => decode_f32(src),
        Type::String => decode_string(src).map(Value::String),
        Type::Hex => decode_hex(src),
        Type::Array => decode_array(src).map(Value::Array),
    }
}

fn decode_i16<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    src.read_i16::<LittleEndian>().map(Value::Int16)
}

fn decode_u16<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    src.read_u16::<LittleEndian>().map(Value::UInt16)
}

fn decode_i32<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    src.read_i32::<LittleEndian>().map(Value::Int32)
}

fn decode_u32<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    src.read_u32::<LittleEndian>().map(Value::UInt32)
}

fn decode_f32<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    src.read_f32::<LittleEndian>().map(Value::Float)
}

fn decode_string<'a>(src: &mut &'a [u8]) -> io::Result<&'a BStr> {
    const NUL: u8 = 0x00;

    let len = src
        .as_bstr()
        .find_byte(NUL)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "string not NUL terminated"))?;

    let (buf, rest) = src.split_at(len);

    *src = &rest[1..];

    Ok(buf.as_bstr())
}

fn decode_hex<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    decode_string(src).map(Value::Hex)
}

fn read_u8(src: &mut &[u8]) -> io::Result<u8> {
    let Some((n, rest)) = src.split_first() else {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    };

    *src = rest;

    Ok(*n)
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
        if let Value::Array(Array::UInt8(values)) = decode_value(&mut src, Type::Array)? {
            let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
            assert_eq!(actual, [0]);
        } else {
            panic!();
        }

        Ok(())
    }
}
