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
    InvalidLength(lexical_core::Error),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidKind(e) => Some(e),
            Self::InvalidLength(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidKind(_) => write!(f, "invalid kind"),
            Self::InvalidLength(_) => write!(f, "invalid length"),
        }
    }
}

pub(crate) fn parse_op(src: &mut &[u8]) -> Result<Op, ParseError> {
    let len = parse_len(src)?;
    let kind = parse_kind(src).map_err(ParseError::InvalidKind)?;
    Ok(Op::new(kind, len))
}

fn parse_len(src: &mut &[u8]) -> Result<usize, ParseError> {
    let (len, i) = lexical_core::parse_partial(src).map_err(ParseError::InvalidLength)?;
    *src = &src[i..];
    Ok(len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_op() {
        use crate::record::cigar::op::Kind;

        fn t(mut src: &[u8], expected: Op) {
            assert_eq!(parse_op(&mut src), Ok(expected));
        }

        t(b"8M", Op::new(Kind::Match, 8));
        t(b"13N", Op::new(Kind::Skip, 13));
        t(b"144S", Op::new(Kind::SoftClip, 144));

        let data = [];
        let mut src = &data[..];
        assert!(matches!(
            parse_op(&mut src),
            Err(ParseError::InvalidLength(_))
        ));

        let data = b"8Z";
        let mut src = &data[..];
        assert!(matches!(
            parse_op(&mut src),
            Err(ParseError::InvalidKind(_))
        ));
    }
}
