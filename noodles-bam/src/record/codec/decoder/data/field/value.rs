mod array;

use std::{error, fmt, string};

use bstr::BString;
use memchr::memchr;
use noodles_sam::alignment::{record::data::field::Type, record_buf::data::field::Value};

use self::array::get_array;

/// An error when a raw BAM record data field value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The string is not NUL terminated.
    StringNotNulTerminated,
    /// The string is invalid.
    InvalidString(string::FromUtf8Error),
    /// The array subtype is invalid.
    InvalidArray(array::DecodeError),
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidString(e) => Some(e),
            Self::InvalidArray(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::StringNotNulTerminated => write!(f, "string is not NUL terminated"),
            Self::InvalidString(_) => write!(f, "invalid string"),
            Self::InvalidArray(_) => write!(f, "invalid array"),
        }
    }
}

pub fn read_value(src: &mut &[u8], ty: Type) -> Result<Value, DecodeError> {
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
        Type::Array => get_array(src).map_err(DecodeError::InvalidArray),
    }
}

fn read_u8(src: &mut &[u8]) -> Result<u8, DecodeError> {
    let (n, rest) = src.split_first().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(*n)
}

fn read_u16_le(src: &mut &[u8]) -> Result<u16, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(u16::from_le_bytes(*buf))
}

fn read_u32_le(src: &mut &[u8]) -> Result<u32, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(u32::from_le_bytes(*buf))
}

fn read_f32_le(src: &mut &[u8]) -> Result<f32, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(f32::from_le_bytes(*buf))
}

fn read_string(src: &mut &[u8]) -> Result<BString, DecodeError> {
    const NUL: u8 = 0x00;

    let i = memchr(NUL, src).ok_or(DecodeError::StringNotNulTerminated)?;
    let (buf, rest) = src.split_at(i);
    *src = &rest[1..];
    Ok(buf.into())
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record_buf::data::field::value::Array;

    use super::*;

    #[test]
    fn test_read_value() -> Result<(), Box<dyn std::error::Error>> {
        fn t(mut src: &[u8], ty: Type, expected: Value) -> Result<(), DecodeError> {
            let actual = read_value(&mut src, ty)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(b"n", Type::Character, Value::Character(b'n'))?;
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
            Value::from("ndls"),
        )?;
        t(
            &[b'C', b'A', b'F', b'E', 0x00],
            Type::Hex,
            Value::Hex(b"CAFE".into()),
        )?;

        t(
            &[b'c', 0x01, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::Array(Array::Int8(vec![0])),
        )?;
        t(
            &[b'C', 0x01, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::Array(Array::UInt8(vec![0])),
        )?;
        t(
            &[b's', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::Array(Array::Int16(vec![0])),
        )?;
        t(
            &[b'S', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::Array(Array::UInt16(vec![0])),
        )?;
        t(
            &[b'i', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::Array(Array::Int32(vec![0])),
        )?;
        t(
            &[b'I', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::Array(Array::UInt32(vec![0])),
        )?;
        t(
            &[b'f', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::Array(Array::Float(vec![0.0])),
        )?;

        Ok(())
    }
}
