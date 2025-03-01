use std::{error, fmt};

use noodles_sam::alignment::record::Flags;

/// An error when raw BAM record flags fail to parse.
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

pub(super) fn read_flags(src: &mut &[u8]) -> Result<Flags, DecodeError> {
    read_u16_le(src).map(Flags::from)
}

fn read_u16_le(src: &mut &[u8]) -> Result<u16, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(u16::from_le_bytes(*buf))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_flags() {
        let mut src = &[0x00, 0x00][..];
        assert_eq!(read_flags(&mut src), Ok(Flags::empty()));

        let mut src = &[0x04, 0x00][..];
        assert_eq!(read_flags(&mut src), Ok(Flags::UNMAPPED));

        let mut src = &[][..];
        assert_eq!(read_flags(&mut src), Err(DecodeError::UnexpectedEof));

        let mut src = &[0x00][..];
        assert_eq!(read_flags(&mut src), Err(DecodeError::UnexpectedEof));
    }
}
