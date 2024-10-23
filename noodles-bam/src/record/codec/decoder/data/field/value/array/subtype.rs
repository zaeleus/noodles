use std::{error, fmt};

use noodles_sam::alignment::record::data::field::value::array::Subtype;

/// An error when a raw BAM record data field value subtype fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The subtype is invalid.
    Invalid { actual: u8 },
}

impl error::Error for DecodeError {}

impl fmt::Display for DecodeError {
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

pub fn read_subtype(src: &mut &[u8]) -> Result<Subtype, DecodeError> {
    let (n, rest) = src.split_first().ok_or(DecodeError::UnexpectedEof)?;

    *src = rest;

    match *n {
        b'c' => Ok(Subtype::Int8),
        b'C' => Ok(Subtype::UInt8),
        b's' => Ok(Subtype::Int16),
        b'S' => Ok(Subtype::UInt16),
        b'i' => Ok(Subtype::Int32),
        b'I' => Ok(Subtype::UInt32),
        b'f' => Ok(Subtype::Float),
        _ => Err(DecodeError::Invalid { actual: *n }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_subtype() -> Result<(), DecodeError> {
        fn t(mut src: &[u8], expected: Subtype) -> Result<(), DecodeError> {
            assert_eq!(read_subtype(&mut src)?, expected);
            Ok(())
        }

        t(b"c", Subtype::Int8)?;
        t(b"C", Subtype::UInt8)?;
        t(b"s", Subtype::Int16)?;
        t(b"S", Subtype::UInt16)?;
        t(b"i", Subtype::Int32)?;
        t(b"I", Subtype::UInt32)?;
        t(b"f", Subtype::Float)?;

        let data = b"";
        let mut src = &data[..];
        assert_eq!(read_subtype(&mut src), Err(DecodeError::UnexpectedEof));

        let data = b"n";
        let mut src = &data[..];
        assert_eq!(
            read_subtype(&mut src),
            Err(DecodeError::Invalid { actual: b'n' })
        );

        Ok(())
    }
}
