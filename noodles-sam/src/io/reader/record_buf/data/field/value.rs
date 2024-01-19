mod array;

use std::{error, fmt, str};

use self::array::parse_array;
use crate::alignment::{record::data::field::Type, record_buf::data::field::Value};

/// An error when a raw SAM record data field value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The type is invalid.
    InvalidType { actual: Type },
    /// The character is invalid.
    InvalidCharacter,
    /// The integer is invalid.
    InvalidInteger(lexical_core::Error),
    /// The integer value is invalid.
    InvalidIntegerValue,
    /// The float is invalid.
    InvalidFloat(lexical_core::Error),
    /// The string is invalid.
    InvalidString,
    /// The hex is invalid.
    InvalidHex,
    /// The array is invalid.
    InvalidArray(array::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidInteger(e) => Some(e),
            Self::InvalidFloat(e) => Some(e),
            Self::InvalidArray(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::InvalidType { actual } => write!(
                f,
                "invalid type: expected {{Character, Int32, Float, String, Hex, Array}}, got {:?}",
                actual
            ),
            Self::InvalidCharacter => write!(f, "invalid character"),
            Self::InvalidInteger(_) => write!(f, "invalid integer"),
            Self::InvalidIntegerValue => write!(f, "invalid integer value"),
            Self::InvalidFloat(_) => write!(f, "invalid float"),
            Self::InvalidString => write!(f, "invalid string"),
            Self::InvalidHex => write!(f, "invalid hex"),
            Self::InvalidArray(_) => write!(f, "invalid array"),
        }
    }
}

pub(super) fn parse_value(src: &mut &[u8], ty: Type) -> Result<Value, ParseError> {
    match ty {
        Type::Character => parse_char(src),
        Type::Int32 => parse_int(src),
        Type::Float => parse_float(src),
        Type::String => parse_string(src),
        Type::Hex => parse_hex(src),
        Type::Array => parse_array(src)
            .map(Value::Array)
            .map_err(ParseError::InvalidArray),
        _ => Err(ParseError::InvalidType { actual: ty }),
    }
}

fn parse_char(src: &[u8]) -> Result<Value, ParseError> {
    let (n, rest) = src.split_first().ok_or(ParseError::UnexpectedEof)?;

    if rest.is_empty() {
        Ok(Value::Character(*n))
    } else {
        Err(ParseError::InvalidCharacter)
    }
}

fn parse_int(src: &[u8]) -> Result<Value, ParseError> {
    lexical_core::parse::<i64>(src)
        .map_err(ParseError::InvalidInteger)
        .and_then(Value::try_from)
}

fn parse_float(src: &[u8]) -> Result<Value, ParseError> {
    lexical_core::parse(src)
        .map(Value::Float)
        .map_err(ParseError::InvalidFloat)
}

fn parse_string(src: &[u8]) -> Result<Value, ParseError> {
    if src.iter().all(|n| matches!(n, b' '..=b'~')) {
        str::from_utf8(src)
            .map(|s| Value::String(s.into()))
            .map_err(|_| ParseError::InvalidString)
    } else {
        Err(ParseError::InvalidString)
    }
}

fn parse_hex(src: &[u8]) -> Result<Value, ParseError> {
    fn is_even(n: usize) -> bool {
        n % 2 == 0
    }

    fn is_upper_ascii_hexdigit(n: u8) -> bool {
        matches!(n, b'0'..=b'9' | b'A'..=b'F')
    }

    if is_even(src.len()) && src.iter().copied().all(is_upper_ascii_hexdigit) {
        Ok(Value::Hex(src.into()))
    } else {
        Err(ParseError::InvalidHex)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_value() -> Result<(), Box<dyn std::error::Error>> {
        use crate::alignment::record_buf::data::field::value::Array;

        fn t(mut src: &[u8], ty: Type, expected: Value) {
            assert_eq!(parse_value(&mut src, ty), Ok(expected));
        }

        t(b"n", Type::Character, Value::Character(b'n'));
        assert_eq!(
            parse_value(&mut &b""[..], Type::Character),
            Err(ParseError::UnexpectedEof)
        );
        assert_eq!(
            parse_value(&mut &b"ndls"[..], Type::Character),
            Err(ParseError::InvalidCharacter)
        );

        t(b"0", Type::Int32, Value::UInt8(0));
        assert!(matches!(
            parse_value(&mut &b""[..], Type::Int32),
            Err(ParseError::InvalidInteger(_))
        ));
        assert!(matches!(
            parse_value(&mut &b"ndls"[..], Type::Int32),
            Err(ParseError::InvalidInteger(_))
        ));

        t(b"0", Type::Float, Value::Float(0.0));
        assert!(matches!(
            parse_value(&mut &b""[..], Type::Float),
            Err(ParseError::InvalidFloat(_))
        ));
        assert!(matches!(
            parse_value(&mut &b"ndls"[..], Type::Float),
            Err(ParseError::InvalidFloat(_))
        ));

        t(b"", Type::String, Value::from(""));
        t(b" ", Type::String, Value::from(" "));
        t(b"ndls", Type::String, Value::from("ndls"));
        assert_eq!(
            parse_value(&mut &[0xf0, 0x9f, 0x8d, 0x9c][..], Type::String),
            Err(ParseError::InvalidString)
        );

        t(b"CAFE", Type::Hex, Value::Hex(b"CAFE".into()));
        assert_eq!(
            parse_value(&mut &b"cafe"[..], Type::Hex),
            Err(ParseError::InvalidHex)
        );
        assert_eq!(
            parse_value(&mut &b"CAFE0"[..], Type::Hex),
            Err(ParseError::InvalidHex)
        );
        assert_eq!(
            parse_value(&mut &b"NDLS"[..], Type::Hex),
            Err(ParseError::InvalidHex)
        );

        t(b"C,0", Type::Array, Value::Array(Array::UInt8(vec![0])));

        Ok(())
    }
}
