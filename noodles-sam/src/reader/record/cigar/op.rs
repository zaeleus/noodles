mod kind;

use std::{error, fmt};

use self::kind::parse_kind;
use crate::record::cigar::Op;

/// An error when a raw SAM record CIGAR op fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The kind is invalid.
    InvalidKind(kind::ParseError),
    /// The length is invalid.
    InvalidLength,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidKind(e) => Some(e),
            Self::InvalidLength => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidKind(_) => write!(f, "invalid kind"),
            Self::InvalidLength => write!(f, "invalid length"),
        }
    }
}

pub(super) fn parse_op(src: &mut &[u8]) -> Result<Op, ParseError> {
    let len = parse_len(src)?;
    let kind = parse_kind(src).map_err(ParseError::InvalidKind)?;
    Ok(Op::new(kind, len))
}

fn parse_len(src: &mut &[u8]) -> Result<usize, ParseError> {
    let (len, i) = lexical_core::parse_partial(src).map_err(|_| ParseError::InvalidLength)?;
    *src = &src[i..];
    Ok(len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_len() {
        assert_eq!(parse_len(&mut &[][..]), Err(ParseError::InvalidLength));
    }
}
