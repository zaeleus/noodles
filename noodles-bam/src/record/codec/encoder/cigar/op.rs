mod kind;

use std::{error, fmt};

use noodles_sam::alignment::record::cigar::Op;

use self::kind::encode_kind;

const MAX_LENGTH: usize = (1 << 28) - 1;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EncodeError {
    InvalidLength(usize),
}

impl error::Error for EncodeError {}

impl fmt::Display for EncodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLength(actual) => {
                write!(f, "invalid length: expected <= {MAX_LENGTH}, got {actual}")
            }
        }
    }
}

pub(super) fn encode_op(op: Op) -> Result<u32, EncodeError> {
    if op.len() <= MAX_LENGTH {
        // SAFETY: `op.len() <= MAX_LENGTH <= u32::MAX`.
        let len = op.len() as u32;
        let k = encode_kind(op.kind());
        Ok((len << 4) | k)
    } else {
        Err(EncodeError::InvalidLength(op.len()))
    }
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record::cigar::op::Kind;

    use super::*;

    #[test]
    fn test_encode_op() {
        assert_eq!(encode_op(Op::new(Kind::Match, 1)), Ok(0x10));

        assert_eq!(
            encode_op(Op::new(Kind::Match, 1 << 28)),
            Err(EncodeError::InvalidLength(1 << 28))
        );
    }
}
