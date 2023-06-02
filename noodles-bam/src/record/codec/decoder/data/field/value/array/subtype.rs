use std::{error, fmt, mem};

use bytes::Buf;
use noodles_sam::record::data::field::value::array::Subtype;

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

pub fn get_subtype<B>(src: &mut B) -> Result<Subtype, DecodeError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u8>() {
        return Err(DecodeError::UnexpectedEof);
    }

    match src.get_u8() {
        b'c' => Ok(Subtype::Int8),
        b'C' => Ok(Subtype::UInt8),
        b's' => Ok(Subtype::Int16),
        b'S' => Ok(Subtype::UInt16),
        b'i' => Ok(Subtype::Int32),
        b'I' => Ok(Subtype::UInt32),
        b'f' => Ok(Subtype::Float),
        n => Err(DecodeError::Invalid { actual: n }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_subtype() -> Result<(), DecodeError> {
        fn t(mut src: &[u8], expected: Subtype) -> Result<(), DecodeError> {
            assert_eq!(get_subtype(&mut src)?, expected);
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
        assert_eq!(get_subtype(&mut src), Err(DecodeError::UnexpectedEof));

        let data = b"n";
        let mut src = &data[..];
        assert_eq!(
            get_subtype(&mut src),
            Err(DecodeError::Invalid { actual: b'n' })
        );

        Ok(())
    }
}
