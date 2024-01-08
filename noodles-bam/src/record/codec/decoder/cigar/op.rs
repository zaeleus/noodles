mod kind;

use std::{error, fmt, num};

use noodles_sam::alignment::record::cigar::Op;

use self::kind::decode_kind;

/// An error when a raw BAM record CIGAR op fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// The kind is invalid.
    InvalidKind(kind::DecodeError),
    /// The length is invalid.
    InvalidLength(num::TryFromIntError),
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidKind(e) => Some(e),
            Self::InvalidLength(e) => Some(e),
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidKind(_) => write!(f, "invalid kind"),
            Self::InvalidLength(_) => write!(f, "invalid length"),
        }
    }
}

pub(crate) fn decode_op(n: u32) -> Result<Op, DecodeError> {
    let kind = decode_kind(n).map_err(DecodeError::InvalidKind)?;
    let len = usize::try_from(n >> 4).map_err(DecodeError::InvalidLength)?;
    Ok(Op::new(kind, len))
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record::cigar::op::Kind;

    use super::*;

    #[test]
    fn test_decode_op() {
        assert_eq!(decode_op(0x10), Ok(Op::new(Kind::Match, 1)));
        assert!(matches!(decode_op(0x19), Err(DecodeError::InvalidKind(_))));
    }
}
