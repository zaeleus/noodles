//! SAM record data field value.

mod array;
pub mod base_modifications;
mod integer;

pub use self::base_modifications::BaseModifications;

use std::io;

use bstr::{BStr, ByteSlice};

use self::{array::parse_array, integer::parse_integer_value};
use super::Type;
use crate::alignment::record::data::field::Value;

pub(super) fn parse_value<'a>(src: &mut &'a [u8], ty: Type) -> io::Result<Value<'a>> {
    match ty {
        Type::Character => parse_character_value(src),
        Type::Integer => parse_integer_value(src),
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

fn parse_float_value<'a>(src: &mut &'a [u8]) -> io::Result<Value<'a>> {
    let (n, i) = lexical_core::parse_partial(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    *src = &src[i..];

    Ok(Value::Float(n))
}

fn parse_string<'a>(src: &mut &'a [u8]) -> &'a BStr {
    const DELIMITER: u8 = b'\t';

    let i = src.as_bstr().find_byte(DELIMITER).unwrap_or(src.len());
    let (buf, rest) = src.split_at(i);
    *src = rest;
    buf.as_bstr()
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
        let mut src = &b"n"[..];
        assert!(matches!(
            parse_value(&mut src, Type::Character)?,
            Value::Character(b'n')
        ));

        let mut src = &b"0"[..];
        assert!(matches!(
            parse_value(&mut src, Type::Integer)?,
            Value::Int32(0)
        ));

        let mut src = &b"0.0"[..];
        assert!(matches!(
            parse_value(&mut src, Type::Float)?,
            Value::Float(n) if n == 0.0
        ));

        let mut src = &b"ndls"[..];
        assert!(matches!(
            parse_value(&mut src, Type::String)?,
            Value::String(s) if s == "ndls"
        ));

        let mut src = &b"CAFE"[..];
        assert!(matches!(
            parse_value(&mut src, Type::Hex)?,
            Value::Hex(s) if s == "CAFE"
        ));

        let mut src = &b"C,0"[..];
        assert!(matches!(
            parse_value(&mut src, Type::Array)?,
            Value::Array(_)
        ));

        Ok(())
    }
}
