use std::{error, fmt};

use noodles_sam::alignment::record::data::field::Tag;

/// An error when a raw BAM record data field tag fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// Unexpected EOF.
    UnexpectedEof,
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::UnexpectedEof => None,
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
        }
    }
}

pub fn read_tag(src: &mut &[u8]) -> Result<Tag, DecodeError> {
    let ([b0, b1], rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(Tag::new(*b0, *b1))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_tag() {
        let data = [b'N', b'H'];
        let mut reader = &data[..];
        assert_eq!(read_tag(&mut reader), Ok(Tag::ALIGNMENT_HIT_COUNT));
    }
}
