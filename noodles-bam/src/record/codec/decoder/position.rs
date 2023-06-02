use std::{error, fmt, mem, num};

use bytes::Buf;
use noodles_core::Position;

/// An error when raw BAM record flags fail to parse.
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

pub(crate) fn get_position<B>(src: &mut B) -> Result<Option<Position>, DecodeError>
where
    B: Buf,
{
    const MISSING: i32 = -1;

    if src.remaining() < mem::size_of::<i32>() {
        return Err(DecodeError::UnexpectedEof);
    }

    match src.get_i32_le() {
        MISSING => Ok(None),
        n => usize::try_from(n)
            .map(|m| m + 1)
            .map_err(DecodeError::Invalid)
            .map(Position::new),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_position() {
        let data = (-1i32).to_le_bytes();
        let mut src = &data[..];
        assert_eq!(get_position(&mut src), Ok(None));

        let data = 0i32.to_le_bytes();
        let mut src = &data[..];
        assert_eq!(get_position(&mut src), Ok(Some(Position::MIN)));

        let mut src = &[][..];
        assert_eq!(get_position(&mut src), Err(DecodeError::UnexpectedEof));

        let data = (-2i32).to_le_bytes();
        let mut src = &data[..];
        assert!(matches!(
            get_position(&mut src),
            Err(DecodeError::Invalid(_))
        ));
    }

    #[cfg(not(target_pointer_width = "16"))]
    #[test]
    fn test_get_position_with_max_position() -> Result<(), num::TryFromIntError> {
        let data = i32::MAX.to_le_bytes();
        let mut src = &data[..];
        let expected = Position::try_from(1 << 31)?;
        assert_eq!(get_position(&mut src), Ok(Some(expected)));
        Ok(())
    }
}
