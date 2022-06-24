mod subtype;
mod ty;

pub(super) use self::ty::parse_type;

use std::{io, str};

use self::subtype::parse_subtype;
use crate::record::data::field::{
    value::{Character, Subtype, Type},
    Value,
};

pub(super) fn parse_value(src: &mut &[u8], ty: Type) -> io::Result<Value> {
    match ty {
        Type::Character => parse_char_value(src),
        Type::Int32 => parse_int_value(src),
        Type::Float => parse_float_value(src),
        Type::String => parse_string_value(src),
        Type::Hex => parse_hex_value(src),
        Type::Array => parse_array_value(src),
        _ => Err(io::Error::from(io::ErrorKind::InvalidData)),
    }
}

fn parse_char_value(src: &[u8]) -> io::Result<Value> {
    let (n, rest) = src
        .split_first()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    if rest.is_empty() {
        Character::try_from(*n)
            .map(Value::Character)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    } else {
        Err(io::Error::from(io::ErrorKind::InvalidData))
    }
}

fn parse_int_value(src: &[u8]) -> io::Result<Value> {
    lexical_core::parse::<i32>(src)
        .map(Value::from)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn parse_float_value(src: &[u8]) -> io::Result<Value> {
    lexical_core::parse(src)
        .map(Value::Float)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn parse_string_value(src: &[u8]) -> io::Result<Value> {
    if src.iter().all(|n| matches!(n, b' '..=b'~')) {
        str::from_utf8(src)
            .map(|s| Value::String(s.into()))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    } else {
        Err(io::Error::from(io::ErrorKind::InvalidData))
    }
}

fn parse_hex_value(src: &[u8]) -> io::Result<Value> {
    if src.len() % 2 == 0 && src.iter().all(|n| n.is_ascii_hexdigit()) {
        str::from_utf8(src)
            .map(|s| Value::Hex(s.into()))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    } else {
        Err(io::Error::from(io::ErrorKind::InvalidData))
    }
}

fn parse_array_value(src: &mut &[u8]) -> io::Result<Value> {
    const DELIMITER: u8 = b',';

    fn consume_delimiter(src: &mut &[u8]) -> io::Result<()> {
        let (n, rest) = src
            .split_first()
            .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

        *src = rest;

        if *n == DELIMITER {
            Ok(())
        } else {
            Err(io::Error::from(io::ErrorKind::InvalidData))
        }
    }

    let subtype = parse_subtype(src)?;

    match subtype {
        Subtype::Int8 => {
            let mut values = Vec::new();

            while !src.is_empty() {
                consume_delimiter(src)?;
                let (value, i) = lexical_core::parse_partial(src)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                *src = &src[i..];
                values.push(value);
            }

            Ok(Value::Int8Array(values))
        }
        Subtype::UInt8 => {
            let mut values = Vec::new();

            while !src.is_empty() {
                consume_delimiter(src)?;
                let (value, i) = lexical_core::parse_partial(src)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                *src = &src[i..];
                values.push(value);
            }

            Ok(Value::UInt8Array(values))
        }
        Subtype::Int16 => {
            let mut values = Vec::new();

            while !src.is_empty() {
                consume_delimiter(src)?;
                let (value, i) = lexical_core::parse_partial(src)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                *src = &src[i..];
                values.push(value);
            }

            Ok(Value::Int16Array(values))
        }
        Subtype::UInt16 => {
            let mut values = Vec::new();

            while !src.is_empty() {
                consume_delimiter(src)?;
                let (value, i) = lexical_core::parse_partial(src)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                *src = &src[i..];
                values.push(value);
            }

            Ok(Value::UInt16Array(values))
        }
        Subtype::Int32 => {
            let mut values = Vec::new();

            while !src.is_empty() {
                consume_delimiter(src)?;
                let (value, i) = lexical_core::parse_partial(src)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                *src = &src[i..];
                values.push(value);
            }

            Ok(Value::Int32Array(values))
        }
        Subtype::UInt32 => {
            let mut values = Vec::new();

            while !src.is_empty() {
                consume_delimiter(src)?;
                let (value, i) = lexical_core::parse_partial(src)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                *src = &src[i..];
                values.push(value);
            }

            Ok(Value::UInt32Array(values))
        }
        Subtype::Float => {
            let mut values = Vec::new();

            while !src.is_empty() {
                consume_delimiter(src)?;
                let (value, i) = lexical_core::parse_partial(src)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                *src = &src[i..];
                values.push(value);
            }

            Ok(Value::FloatArray(values))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_value() -> Result<(), Box<dyn std::error::Error>> {
        fn t(mut src: &[u8], ty: Type, expected: Value) -> io::Result<()> {
            let actual = parse_value(&mut src, ty)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(
            b"n",
            Type::Character,
            Value::Character(Character::try_from('n')?),
        )?;
        t(b"0", Type::Int32, Value::UInt8(0))?;
        t(b"0", Type::Float, Value::Float(0.0))?;
        t(b"ndls", Type::String, Value::String(String::from("ndls")))?;
        t(b"CAFE", Type::Hex, Value::Hex(String::from("CAFE")))?;

        t(b"c", Type::Array, Value::Int8Array(vec![]))?;
        t(b"c,0", Type::Array, Value::Int8Array(vec![0]))?;
        t(b"c,0,0", Type::Array, Value::Int8Array(vec![0, 0]))?;

        t(b"C", Type::Array, Value::UInt8Array(vec![]))?;
        t(b"C,0", Type::Array, Value::UInt8Array(vec![0]))?;
        t(b"C,0,0", Type::Array, Value::UInt8Array(vec![0, 0]))?;

        t(b"s", Type::Array, Value::Int16Array(vec![]))?;
        t(b"s,0", Type::Array, Value::Int16Array(vec![0]))?;
        t(b"s,0,0", Type::Array, Value::Int16Array(vec![0, 0]))?;

        t(b"S", Type::Array, Value::UInt16Array(vec![]))?;
        t(b"S,0", Type::Array, Value::UInt16Array(vec![0]))?;
        t(b"S,0,0", Type::Array, Value::UInt16Array(vec![0, 0]))?;

        t(b"i", Type::Array, Value::Int32Array(vec![]))?;
        t(b"i,0", Type::Array, Value::Int32Array(vec![0]))?;
        t(b"i,0,0", Type::Array, Value::Int32Array(vec![0, 0]))?;

        t(b"I", Type::Array, Value::UInt32Array(vec![]))?;
        t(b"I,0", Type::Array, Value::UInt32Array(vec![0]))?;
        t(b"I,0,0", Type::Array, Value::UInt32Array(vec![0, 0]))?;

        t(b"f", Type::Array, Value::FloatArray(vec![]))?;
        t(b"f,0", Type::Array, Value::FloatArray(vec![0.0]))?;
        t(b"f,0,0", Type::Array, Value::FloatArray(vec![0.0, 0.0]))?;

        Ok(())
    }
}
