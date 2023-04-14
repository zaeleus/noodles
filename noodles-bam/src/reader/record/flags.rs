use std::{error, fmt, mem};

use bytes::Buf;
use noodles_sam::record::Flags;

/// An error when raw BAM record flags fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Unexpected EOF.
    UnexpectedEof,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
        }
    }
}

pub(crate) fn get_flags<B>(src: &mut B) -> Result<Flags, ParseError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u16>() {
        return Err(ParseError::UnexpectedEof);
    }

    Ok(Flags::from(src.get_u16_le()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_flags() {
        let mut src = &[0x00, 0x00][..];
        assert_eq!(get_flags(&mut src), Ok(Flags::empty()));

        let mut src = &[0x04, 0x00][..];
        assert_eq!(get_flags(&mut src), Ok(Flags::UNMAPPED));

        let mut src = &[][..];
        assert_eq!(get_flags(&mut src), Err(ParseError::UnexpectedEof));

        let mut src = &[0x00][..];
        assert_eq!(get_flags(&mut src), Err(ParseError::UnexpectedEof));
    }
}
