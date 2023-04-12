mod array;

use std::{io, str};

use self::array::parse_array;
use crate::record::data::field::{
    value::{Character, Hex},
    Type, Value,
};

pub(crate) fn parse_value(src: &mut &[u8], ty: Type) -> io::Result<Value> {
    match ty {
        Type::Character => parse_char(src),
        Type::Int32 => parse_int(src),
        Type::Float => parse_float(src),
        Type::String => parse_string(src),
        Type::Hex => parse_hex(src),
        Type::Array => parse_array(src)
            .map(Value::Array)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
        _ => Err(io::Error::from(io::ErrorKind::InvalidData)),
    }
}

fn parse_char(src: &[u8]) -> io::Result<Value> {
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

fn parse_int(src: &[u8]) -> io::Result<Value> {
    lexical_core::parse::<i32>(src)
        .map(Value::from)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn parse_float(src: &[u8]) -> io::Result<Value> {
    lexical_core::parse(src)
        .map(Value::Float)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn parse_string(src: &[u8]) -> io::Result<Value> {
    if src.iter().all(|n| matches!(n, b' '..=b'~')) {
        str::from_utf8(src)
            .map(|s| Value::String(s.into()))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    } else {
        Err(io::Error::from(io::ErrorKind::InvalidData))
    }
}

fn parse_hex(src: &[u8]) -> io::Result<Value> {
    Hex::try_from(src)
        .map(Value::Hex)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_value() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::data::field::value::Array;

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
        assert!(matches!(
            parse_value(&mut &b""[..], Type::Character),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));
        assert!(matches!(
            parse_value(&mut &b"ndls"[..], Type::Character),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        t(b"0", Type::Int32, Value::UInt8(0))?;
        assert!(matches!(
            parse_value(&mut &b""[..], Type::Int32),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
        assert!(matches!(
            parse_value(&mut &b"ndls"[..], Type::Int32),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        t(b"0", Type::Float, Value::Float(0.0))?;
        assert!(matches!(
            parse_value(&mut &b""[..], Type::Float),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
        assert!(matches!(
            parse_value(&mut &b"ndls"[..], Type::Float),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        t(b"", Type::String, Value::String(String::new()))?;
        t(b" ", Type::String, Value::String(String::from(" ")))?;
        t(b"ndls", Type::String, Value::String(String::from("ndls")))?;
        assert!(matches!(
            parse_value(&mut &[0xf0, 0x9f, 0x8d, 0x9c][..], Type::String),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        t(b"CAFE", Type::Hex, Value::Hex("CAFE".parse()?))?;
        assert!(matches!(
            parse_value(&mut &b"cafe"[..], Type::Hex),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
        assert!(matches!(
            parse_value(&mut &b"CAFE0"[..], Type::Hex),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
        assert!(matches!(
            parse_value(&mut &b"NDLS"[..], Type::Hex),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        t(b"C,0", Type::Array, Value::Array(Array::UInt8(vec![0])))?;

        Ok(())
    }
}
