use std::{error, fmt};

/// An error when a raw BAM record template length fails to parse.
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

pub(super) fn read_template_length(src: &mut &[u8]) -> Result<i32, DecodeError> {
    read_i32_le(src)
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
    fn test_read_template_length() {
        let mut src = &144i32.to_le_bytes()[..];
        assert_eq!(read_template_length(&mut src), Ok(144));

        let mut src = &[][..];
        assert_eq!(
            read_template_length(&mut src),
            Err(DecodeError::UnexpectedEof),
        );
    }
}
