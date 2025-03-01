use std::{error, fmt};

/// An error when a raw BAM record reference sequence ID fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The input is invalid.
    Invalid,
}

impl error::Error for DecodeError {}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::Invalid => write!(f, "invalid input"),
        }
    }
}

pub(super) fn read_reference_sequence_id(src: &mut &[u8]) -> Result<Option<usize>, DecodeError> {
    const UNMAPPED: i32 = -1;

    match read_i32_le(src)? {
        UNMAPPED => Ok(None),
        n => usize::try_from(n)
            .map(Some)
            .map_err(|_| DecodeError::Invalid),
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
    fn test_read_reference_sequence_id() {
        let data = (-1i32).to_le_bytes();
        let mut src = &data[..];
        assert_eq!(read_reference_sequence_id(&mut src), Ok(None));

        let data = 0i32.to_le_bytes();
        let mut src = &data[..];
        assert_eq!(read_reference_sequence_id(&mut src), Ok(Some(0)));

        let data = [];
        let mut src = &data[..];
        assert_eq!(
            read_reference_sequence_id(&mut src),
            Err(DecodeError::UnexpectedEof)
        );

        let data = (-2i32).to_le_bytes();
        let mut src = &data[..];
        assert_eq!(
            read_reference_sequence_id(&mut src),
            Err(DecodeError::Invalid)
        );
    }
}
