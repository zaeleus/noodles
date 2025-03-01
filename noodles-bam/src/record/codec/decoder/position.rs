use std::{error, fmt, num};

use noodles_core::Position;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The input is invalid.
    Invalid(num::TryFromIntError),
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::UnexpectedEof => None,
            Self::Invalid(e) => Some(e),
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::Invalid(_) => write!(f, "invalid input"),
        }
    }
}

pub(super) fn read_position(src: &mut &[u8]) -> Result<Option<Position>, DecodeError> {
    const MISSING: i32 = -1;

    match read_i32_le(src)? {
        MISSING => Ok(None),
        n => usize::try_from(n)
            .map(|m| m + 1)
            .map_err(DecodeError::Invalid)
            .map(Position::new),
    }
}

fn read_i32_le(src: &mut &[u8]) -> Result<i32, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(i32::from_le_bytes(*buf))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_position() {
        let data = (-1i32).to_le_bytes();
        let mut src = &data[..];
        assert_eq!(read_position(&mut src), Ok(None));

        let data = 0i32.to_le_bytes();
        let mut src = &data[..];
        assert_eq!(read_position(&mut src), Ok(Some(Position::MIN)));

        let mut src = &[][..];
        assert_eq!(read_position(&mut src), Err(DecodeError::UnexpectedEof));

        let data = (-2i32).to_le_bytes();
        let mut src = &data[..];
        assert!(matches!(
            read_position(&mut src),
            Err(DecodeError::Invalid(_))
        ));
    }

    #[cfg(not(target_pointer_width = "16"))]
    #[test]
    fn test_read_position_with_max_position() -> Result<(), num::TryFromIntError> {
        let data = i32::MAX.to_le_bytes();
        let mut src = &data[..];
        let expected = Position::try_from(1 << 31)?;
        assert_eq!(read_position(&mut src), Ok(Some(expected)));
        Ok(())
    }
}
