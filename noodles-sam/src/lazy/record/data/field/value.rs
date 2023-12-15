mod array;

use std::io;

use self::array::{parse_array, Array};
use super::Type;

/// A raw SAM record data field value.
#[derive(Debug, PartialEq)]
pub enum Value<'a> {
    /// A character (`A`).
    Character(u8),
    /// A 32-bit integer (`i`).
    Int32(i32),
    /// A single-precision floating-point (`f`).
    Float(f32),
    /// A string (`Z`).
    String(&'a [u8]),
    /// A hex string (`H`).
    Hex(&'a [u8]),
    /// An array (`B`).
    Array(Array<'a>),
}

pub(super) fn parse_value<'a>(src: &mut &'a [u8], ty: Type) -> io::Result<Value<'a>> {
    match ty {
        Type::Character => parse_character_value(src),
        Type::Int32 => parse_int32_value(src),
        Type::Float => parse_float_value(src),
        Type::String => Ok(parse_string_value(src)),
        Type::Hex => Ok(parse_hex_value(src)),
        Type::Array => parse_array(src).map(Value::Array),
    }
}

fn parse_character_value<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    if let Some((b, rest)) = src.split_first() {
        *src = rest;
        Ok(Value::Character(*b))
    } else {
        Err(io::Error::from(io::ErrorKind::UnexpectedEof))
    }
}

fn parse_int32_value<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    let (n, i) = lexical_core::parse_partial(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *src = &src[i..];

    Ok(Value::Int32(n))
}

fn parse_float_value<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    let (n, i) = lexical_core::parse_partial(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *src = &src[i..];

    Ok(Value::Float(n))
}

fn parse_string<'a>(src: &mut &'a [u8]) -> &'a [u8] {
    const DELIMITER: u8 = b'\t';

    let i = src
        .iter()
        .position(|&b| b == DELIMITER)
        .unwrap_or(src.len());

    let (buf, rest) = src.split_at(i);

    *src = rest;

    buf
}

fn parse_string_value<'a>(src: &mut &'a [u8]) -> Value<'a> {
    Value::String(parse_string(src))
}

fn parse_hex_value<'a>(src: &mut &'a [u8]) -> Value<'a> {
    Value::Hex(parse_string(src))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_value() -> io::Result<()> {
        fn t(mut data: &[u8], ty: Type, expected: Value<'_>) -> io::Result<()> {
            assert_eq!(parse_value(&mut data, ty)?, expected);
            Ok(())
        }

        t(&b"n"[..], Type::Character, Value::Character(b'n'))?;
        t(&b"0"[..], Type::Int32, Value::Int32(0))?;
        t(&b"0.0"[..], Type::Float, Value::Float(0.0))?;
        t(&b"ndls"[..], Type::String, Value::String(b"ndls"))?;
        t(&b"CAFE"[..], Type::Hex, Value::Hex(b"CAFE"))?;

        t(
            &b"C,0"[..],
            Type::Array,
            Value::Array(Array::UInt8(array::Values::new(b"0"))),
        )?;

        Ok(())
    }
}
