//! Raw BAM record data field value.

mod array;

use std::io::{self, BufRead};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_sam::record::data::field::Type;

use self::array::decode_array;
pub use self::array::Array;

/// A raw BAM record data field value.
#[derive(Debug, PartialEq)]
pub enum Value<'a> {
    /// A character (`A`).
    Character(u8),
    /// An 8-bit integer (`c`).
    Int8(i8),
    /// An 8-bit unsigned integer (`C`).
    UInt8(u8),
    /// A 16-bit integer (`s`).
    Int16(i16),
    /// A 16-bit unsigned integer (`S`).
    UInt16(u16),
    /// A 32-bit integer (`i`).
    Int32(i32),
    /// A 32-bit unsigned integer (`I`).
    UInt32(u32),
    /// A single-precision floating-point (`f`).
    Float(f32),
    /// A string (`Z`).
    String(&'a [u8]),
    /// A hex string (`H`).
    Hex(&'a [u8]),
    /// An array (`B`).
    Array(Array<'a>),
}

impl<'a> Value<'a> {
    /// Returns the type of the value.
    pub fn ty(&self) -> Type {
        match self {
            Self::Character(_) => Type::Character,
            Self::Int8(_) => Type::Int8,
            Self::UInt8(_) => Type::UInt8,
            Self::Int16(_) => Type::Int16,
            Self::UInt16(_) => Type::UInt16,
            Self::Int32(_) => Type::Int32,
            Self::UInt32(_) => Type::UInt32,
            Self::Float(_) => Type::Float,
            Self::String(_) => Type::String,
            Self::Hex(_) => Type::Hex,
            Self::Array(_) => Type::Array,
        }
    }

    /// Returns the value as a 64-bit integer.
    ///
    /// This is a convenience method that converts any integer to an `i64`, which captures the
    /// entire range of all record data field integer values.
    /// ```
    pub fn as_int(&self) -> Option<i64> {
        match *self {
            Self::Int8(n) => Some(i64::from(n)),
            Self::UInt8(n) => Some(i64::from(n)),
            Self::Int16(n) => Some(i64::from(n)),
            Self::UInt16(n) => Some(i64::from(n)),
            Self::Int32(n) => Some(i64::from(n)),
            Self::UInt32(n) => Some(i64::from(n)),
            _ => None,
        }
    }
}

pub(super) fn decode_value<'a>(src: &mut &'a [u8], ty: Type) -> io::Result<Value<'a>> {
    match ty {
        Type::Character => decode_character(src),
        Type::Int8 => decode_i8(src),
        Type::UInt8 => decode_u8(src),
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

fn decode_character<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    src.read_u8().map(Value::Character)
}

fn decode_i8<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    src.read_i8().map(Value::Int8)
}

fn decode_u8<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    src.read_u8().map(Value::UInt8)
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

fn decode_string<'a>(src: &mut &'a [u8]) -> io::Result<&'a [u8]> {
    const NUL: u8 = 0x00;

    let len = src
        .iter()
        .position(|&b| b == NUL)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "string not NUL terminated"))?;

    let buf = &src[..len];

    // +1 for the terminator.
    src.consume(len + 1);

    Ok(buf)
}

fn decode_hex<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    decode_string(src).map(Value::Hex)
}

#[cfg(test)]
mod tests {
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
        assert_eq!(Value::String(b"ndls").ty(), Type::String);
        assert_eq!(Value::Hex(b"CAFE").ty(), Type::Hex);
        assert_eq!(Value::Array(Array::UInt8(&[0])).ty(), Type::Array);
    }

    #[test]
    fn test_decode_value() -> io::Result<()> {
        fn t(mut data: &[u8], ty: Type, expected: Value<'_>) -> io::Result<()> {
            assert_eq!(decode_value(&mut data, ty)?, expected);
            Ok(())
        }

        t(&[b'n'], Type::Character, Value::Character(b'n'))?;
        t(&[0x00], Type::Int8, Value::Int8(0))?;
        t(&[0x00], Type::UInt8, Value::UInt8(0))?;
        t(&[0x00, 0x00], Type::Int16, Value::Int16(0))?;
        t(&[0x00, 0x00], Type::UInt16, Value::UInt16(0))?;
        t(&[0x00, 0x00, 0x00, 0x00], Type::Int32, Value::Int32(0))?;
        t(&[0x00, 0x00, 0x00, 0x00], Type::UInt32, Value::UInt32(0))?;
        t(&[0x00, 0x00, 0x00, 0x00], Type::Float, Value::Float(0.0))?;
        t(
            &[b'n', b'd', b'l', b's', 0x00],
            Type::String,
            Value::String(b"ndls"),
        )?;
        t(
            &[b'C', b'A', b'F', b'E', 0x00],
            Type::Hex,
            Value::Hex(b"CAFE"),
        )?;

        t(
            &[b'C', 0x01, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::Array(Array::UInt8(&[0x00])),
        )?;

        Ok(())
    }
}
