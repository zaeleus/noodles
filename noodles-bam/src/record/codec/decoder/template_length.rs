use std::{error, fmt, mem};

use bytes::Buf;

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

pub(super) fn get_template_length<B>(src: &mut B) -> Result<i32, DecodeError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<i32>() {
        return Err(DecodeError::UnexpectedEof);
    }

    Ok(src.get_i32_le())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_template_length() {
        let mut src = &144i32.to_le_bytes()[..];
        assert_eq!(get_template_length(&mut src), Ok(144));

        let mut src = &[][..];
        assert_eq!(
            get_template_length(&mut src),
            Err(DecodeError::UnexpectedEof),
        );
    }
}
