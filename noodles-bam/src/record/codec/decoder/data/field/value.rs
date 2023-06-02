mod array;

use std::{error, fmt, mem, string};

use bytes::Buf;
use noodles_sam::record::data::field::{
    value::{character, hex, Character},
    Type, Value,
};

use self::array::get_array;

/// An error when a raw BAM record data field value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The character is invalid.
    InvalidCharacter(character::ParseError),
    /// The string is not NUL terminated.
    StringNotNulTerminated,
    /// The string is invalid.
    InvalidString(string::FromUtf8Error),
    /// The hex is invalid.
    InvalidHex(hex::ParseError),
    /// The array subtype is invalid.
    InvalidArray(array::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidCharacter(e) => Some(e),
            Self::InvalidString(e) => Some(e),
            Self::InvalidHex(e) => Some(e),
            Self::InvalidArray(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::InvalidCharacter(_) => write!(f, "invalid character"),
            Self::StringNotNulTerminated => write!(f, "string is not NUL terminated"),
            Self::InvalidString(_) => write!(f, "invalid string"),
            Self::InvalidHex(_) => write!(f, "invalid hex"),
            Self::InvalidArray(_) => write!(f, "invalid array"),
        }
    }
}

pub fn get_value<B>(src: &mut B, ty: Type) -> Result<Value, ParseError>
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
        Type::Array => get_array(src).map_err(ParseError::InvalidArray),
    }
}

fn get_char<B>(src: &mut B) -> Result<Value, ParseError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u8>() {
        return Err(ParseError::UnexpectedEof);
    }

    Character::try_from(src.get_u8())
        .map(Value::Character)
        .map_err(ParseError::InvalidCharacter)
}

fn get_i8<B>(src: &mut B) -> Result<Value, ParseError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<i8>() {
        return Err(ParseError::UnexpectedEof);
    }

    Ok(Value::Int8(src.get_i8()))
}

fn get_u8<B>(src: &mut B) -> Result<Value, ParseError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u8>() {
        return Err(ParseError::UnexpectedEof);
    }

    Ok(Value::UInt8(src.get_u8()))
}

fn get_i16<B>(src: &mut B) -> Result<Value, ParseError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<i16>() {
        return Err(ParseError::UnexpectedEof);
    }

    Ok(Value::Int16(src.get_i16_le()))
}

fn get_u16<B>(src: &mut B) -> Result<Value, ParseError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u16>() {
        return Err(ParseError::UnexpectedEof);
    }

    Ok(Value::UInt16(src.get_u16_le()))
}

fn get_i32<B>(src: &mut B) -> Result<Value, ParseError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<i32>() {
        return Err(ParseError::UnexpectedEof);
    }

    Ok(Value::Int32(src.get_i32_le()))
}

fn get_u32<B>(src: &mut B) -> Result<Value, ParseError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u32>() {
        return Err(ParseError::UnexpectedEof);
    }

    Ok(Value::UInt32(src.get_u32_le()))
}

fn get_f32<B>(src: &mut B) -> Result<Value, ParseError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<f32>() {
        return Err(ParseError::UnexpectedEof);
    }

    Ok(Value::Float(src.get_f32_le()))
}

fn get_string<B>(src: &mut B) -> Result<String, ParseError>
where
    B: Buf,
{
    const NUL: u8 = 0x00;

    let len = src
        .chunk()
        .iter()
        .position(|&b| b == NUL)
        .ok_or(ParseError::StringNotNulTerminated)?;

    let mut buf = vec![0; len];
    src.copy_to_slice(&mut buf);
    src.advance(1); // Discard the NUL terminator.

    String::from_utf8(buf).map_err(ParseError::InvalidString)
}

fn get_hex<B>(src: &mut B) -> Result<Value, ParseError>
where
    B: Buf,
{
    get_string(src)
        .and_then(|s| s.parse().map_err(ParseError::InvalidHex))
        .map(Value::Hex)
}

#[cfg(test)]
mod tests {
    use noodles_sam::record::data::field::value::Array;

    use super::*;

    #[test]
    fn test_get_value() -> Result<(), Box<dyn std::error::Error>> {
        fn t(mut data: &[u8], ty: Type, expected: Value) -> Result<(), ParseError> {
            let actual = get_value(&mut data, ty)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(
            &[b'n'],
            Type::Character,
            Value::Character(Character::try_from('n')?),
        )?;
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
            Value::String(String::from("ndls")),
        )?;
        t(
            &[b'C', b'A', b'F', b'E', 0x00],
            Type::Hex,
            Value::Hex("CAFE".parse()?),
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
