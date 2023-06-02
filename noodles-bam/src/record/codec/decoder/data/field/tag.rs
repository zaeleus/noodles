use std::{error, fmt};

use bytes::Buf;
use noodles_sam::record::data::field::{tag, Tag};

/// An error when a raw BAM record data field tag fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The input is invalid.
    Invalid(tag::ParseError),
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

pub fn get_tag<B>(src: &mut B) -> Result<Tag, DecodeError>
where
    B: Buf,
{
    let mut buf = [0; 2];

    if src.remaining() < buf.len() {
        return Err(DecodeError::UnexpectedEof);
    }

    src.copy_to_slice(&mut buf);

    Tag::try_from(buf).map_err(DecodeError::Invalid)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_tag() {
        let data = [b'N', b'H'];
        let mut reader = &data[..];
        assert_eq!(get_tag(&mut reader), Ok(tag::ALIGNMENT_HIT_COUNT));
    }
}
