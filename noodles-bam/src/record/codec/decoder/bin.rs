use std::{error, fmt, mem};

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

pub(super) fn consume_bin(src: &mut &[u8]) -> Result<(), DecodeError> {
    let (_, rest) = src
        .split_at_checked(mem::size_of::<u16>())
        .ok_or(DecodeError::UnexpectedEof)?;

    *src = rest;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_consume_bin() {
        let mut src = &[0x00, 0x00][..];
        assert!(consume_bin(&mut src).is_ok());

        let mut src = &[][..];
        assert_eq!(consume_bin(&mut src), Err(DecodeError::UnexpectedEof));

        let mut src = &[0x00][..];
        assert_eq!(consume_bin(&mut src), Err(DecodeError::UnexpectedEof));
    }
}
