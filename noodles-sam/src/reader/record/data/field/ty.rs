use std::{error, fmt};

use crate::record::data::field::Type;

/// An error when a raw SAM record data field type fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The type is invalid.
    Invalid { actual: u8 },
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::Invalid { actual } => write!(
                f,
                "invalid type: expected {{A, i, f, Z, H, B}}, got {}",
                char::from(*actual)
            ),
        }
    }
}

pub(crate) fn parse_type(src: &mut &[u8]) -> Result<Type, ParseError> {
    let (n, rest) = src.split_first().ok_or(ParseError::UnexpectedEof)?;

    *src = rest;

    match n {
        b'A' => Ok(Type::Character),
        b'i' => Ok(Type::Int32),
        b'f' => Ok(Type::Float),
        b'Z' => Ok(Type::String),
        b'H' => Ok(Type::Hex),
        b'B' => Ok(Type::Array),
        _ => Err(ParseError::Invalid { actual: *n }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_type() -> Result<(), ParseError> {
        fn t(mut src: &[u8], expected: Type) -> Result<(), ParseError> {
            assert_eq!(parse_type(&mut src)?, expected);
            Ok(())
        }

        t(b"A", Type::Character)?;
        t(b"i", Type::Int32)?;
        t(b"f", Type::Float)?;
        t(b"Z", Type::String)?;
        t(b"H", Type::Hex)?;
        t(b"B", Type::Array)?;

        let data = b"";
        let mut src = &data[..];
        assert_eq!(parse_type(&mut src), Err(ParseError::UnexpectedEof));

        for &n in b"cCsSIn" {
            let data = [n];
            let mut src = &data[..];
            assert_eq!(parse_type(&mut src), Err(ParseError::Invalid { actual: n }));
        }

        Ok(())
    }
}
