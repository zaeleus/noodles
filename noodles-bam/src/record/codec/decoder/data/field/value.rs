mod array;

use std::{error, fmt, mem, string};

use bstr::{BString, ByteSlice};
use bytes::Buf;
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

pub fn get_value<B>(src: &mut B, ty: Type) -> Result<Value, DecodeError>
where
    B: Buf,
{
    match ty {
        Type::Character => get_char(src),
        Type::Int8 => get_i8(src),
        Type::UInt8 => get_u8(src),
        Type::Int16 => get_i16(src),
        Type::UInt16 => get_u16(src),
        Type::Int32 => get_i32(src),
        Type::UInt32 => get_u32(src),
        Type::Float => get_f32(src),
        Type::String => get_string(src).map(Value::String),
        Type::Hex => get_hex(src),
        Type::Array => get_array(src).map_err(DecodeError::InvalidArray),
    }
}

fn get_char<B>(src: &mut B) -> Result<Value, DecodeError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u8>() {
        return Err(DecodeError::UnexpectedEof);
    }

    Ok(Value::Character(src.get_u8()))
}

fn get_i8<B>(src: &mut B) -> Result<Value, DecodeError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<i8>() {
        return Err(DecodeError::UnexpectedEof);
    }

    Ok(Value::Int8(src.get_i8()))
}

fn get_u8<B>(src: &mut B) -> Result<Value, DecodeError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u8>() {
        return Err(DecodeError::UnexpectedEof);
    }

    Ok(Value::UInt8(src.get_u8()))
}

fn get_i16<B>(src: &mut B) -> Result<Value, DecodeError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<i16>() {
        return Err(DecodeError::UnexpectedEof);
    }

    Ok(Value::Int16(src.get_i16_le()))
}

fn get_u16<B>(src: &mut B) -> Result<Value, DecodeError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u16>() {
        return Err(DecodeError::UnexpectedEof);
    }

    Ok(Value::UInt16(src.get_u16_le()))
}

fn get_i32<B>(src: &mut B) -> Result<Value, DecodeError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<i32>() {
        return Err(DecodeError::UnexpectedEof);
    }

    Ok(Value::Int32(src.get_i32_le()))
}

fn get_u32<B>(src: &mut B) -> Result<Value, DecodeError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u32>() {
        return Err(DecodeError::UnexpectedEof);
    }

    Ok(Value::UInt32(src.get_u32_le()))
}

fn get_f32<B>(src: &mut B) -> Result<Value, DecodeError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<f32>() {
        return Err(DecodeError::UnexpectedEof);
    }

    Ok(Value::Float(src.get_f32_le()))
}

fn get_string<B>(src: &mut B) -> Result<BString, DecodeError>
where
    B: Buf,
{
    const NUL: u8 = 0x00;

    let len = src
        .chunk()
        .as_bstr()
        .find_byte(NUL)
        .ok_or(DecodeError::StringNotNulTerminated)?;

    let mut buf = vec![0; len];
    src.copy_to_slice(&mut buf);
    src.advance(1); // Discard the NUL terminator.

    Ok(buf.into())
}

fn get_hex<B>(src: &mut B) -> Result<Value, DecodeError>
where
    B: Buf,
{
    get_string(src).map(Value::Hex)
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record_buf::data::field::value::Array;

    use super::*;

    #[test]
    fn test_get_value() -> Result<(), Box<dyn std::error::Error>> {
        fn t(mut data: &[u8], ty: Type, expected: Value) -> Result<(), DecodeError> {
            let actual = get_value(&mut data, ty)?;
            assert_eq!(actual, expected);
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
