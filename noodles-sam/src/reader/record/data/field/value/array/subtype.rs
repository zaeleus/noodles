use std::{error, fmt};

use crate::record::data::field::value::array::Subtype;

/// An error when a raw BAM record data field array value subtype fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The subtype is invalid.
    Invalid { actual: u8 },
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::Invalid { actual } => write!(
                f,
                "invalid input: expected {{c, C, s, S, i, I, f}}, got {}",
                char::from(*actual)
            ),
        }
    }
}

pub(super) fn parse_subtype(src: &mut &[u8]) -> Result<Subtype, ParseError> {
    let (n, rest) = src.split_first().ok_or(ParseError::UnexpectedEof)?;

    *src = rest;

    match n {
        b'c' => Ok(Subtype::Int8),
        b'C' => Ok(Subtype::UInt8),
        b's' => Ok(Subtype::Int16),
        b'S' => Ok(Subtype::UInt16),
        b'i' => Ok(Subtype::Int32),
        b'I' => Ok(Subtype::UInt32),
        b'f' => Ok(Subtype::Float),
        _ => Err(ParseError::Invalid { actual: *n }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_subtype() {
        fn t(mut src: &[u8], expected: Subtype) {
            assert_eq!(parse_subtype(&mut src), Ok(expected));
        }

        t(b"c", Subtype::Int8);
        t(b"C", Subtype::UInt8);
        t(b"s", Subtype::Int16);
        t(b"S", Subtype::UInt16);
        t(b"i", Subtype::Int32);
        t(b"I", Subtype::UInt32);
        t(b"f", Subtype::Float);

        let mut src = &b""[..];
        assert_eq!(parse_subtype(&mut src), Err(ParseError::UnexpectedEof));

        let mut src = &b"n"[..];
        assert!(matches!(
            parse_subtype(&mut src),
            Err(ParseError::Invalid { actual: b'n' })
        ));
    }
}
