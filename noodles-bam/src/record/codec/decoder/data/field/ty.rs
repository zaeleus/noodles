use std::{error, fmt};

use noodles_sam::alignment::record::data::field::Type;

/// An error when a raw BAM record data field type fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The type is invalid.
    Invalid { actual: u8 },
}

impl error::Error for DecodeError {}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::Invalid { actual } => write!(
                f,
                "invalid type: expected {{A, c, C, s, S, i, I, f, Z, H, B}}, got {}",
                char::from(*actual)
            ),
        }
    }
}

pub fn read_type(src: &mut &[u8]) -> Result<Type, DecodeError> {
    let (n, rest) = src.split_first().ok_or(DecodeError::UnexpectedEof)?;

    *src = rest;

    match *n {
        b'A' => Ok(Type::Character),
        b'c' => Ok(Type::Int8),
        b'C' => Ok(Type::UInt8),
        b's' => Ok(Type::Int16),
        b'S' => Ok(Type::UInt16),
        b'i' => Ok(Type::Int32),
        b'I' => Ok(Type::UInt32),
        b'f' => Ok(Type::Float),
        b'Z' => Ok(Type::String),
        b'H' => Ok(Type::Hex),
        b'B' => Ok(Type::Array),
        _ => Err(DecodeError::Invalid { actual: *n }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_type() -> Result<(), DecodeError> {
        fn t(mut src: &[u8], expected: Type) -> Result<(), DecodeError> {
            assert_eq!(read_type(&mut src)?, expected);
            Ok(())
        }

        t(b"A", Type::Character)?;
        t(b"c", Type::Int8)?;
        t(b"C", Type::UInt8)?;
        t(b"s", Type::Int16)?;
        t(b"S", Type::UInt16)?;
        t(b"i", Type::Int32)?;
        t(b"I", Type::UInt32)?;
        t(b"f", Type::Float)?;
        t(b"Z", Type::String)?;
        t(b"H", Type::Hex)?;
        t(b"B", Type::Array)?;

        let data = b"";
        let mut src = &data[..];
        assert_eq!(read_type(&mut src), Err(DecodeError::UnexpectedEof));

        let data = b"n";
        let mut src = &data[..];
        assert_eq!(
            read_type(&mut src),
            Err(DecodeError::Invalid { actual: b'n' })
        );

        Ok(())
    }
}
