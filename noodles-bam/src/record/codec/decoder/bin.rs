use std::{error, fmt, mem};

use bytes::Buf;

/// An error when raw BAM record bin fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// Unexpected EOF.
    UnexpectedEof,
}

impl error::Error for DecodeError {}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
        }
    }
}

pub(super) fn discard_bin<B>(src: &mut B) -> Result<(), DecodeError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u16>() {
        return Err(DecodeError::UnexpectedEof);
    }

    src.advance(mem::size_of::<u16>());

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_discard_bin() {
        let mut src = &[0x00, 0x00][..];
        assert!(discard_bin(&mut src).is_ok());

        let mut src = &[][..];
        assert_eq!(discard_bin(&mut src), Err(DecodeError::UnexpectedEof));

        let mut src = &[0x00][..];
        assert_eq!(discard_bin(&mut src), Err(DecodeError::UnexpectedEof));
    }
}
