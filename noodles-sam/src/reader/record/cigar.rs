mod op;

use std::{error, fmt, mem};

use self::op::parse_op;
use crate::record::Cigar;

/// An error when a raw SAM record CIGAR fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
    /// An op is invalid.
    InvalidOp(op::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Invalid => None,
            Self::InvalidOp(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => write!(f, "invalid input"),
            Self::InvalidOp(_) => write!(f, "invalid op"),
        }
    }
}

pub(crate) fn parse_cigar(mut src: &[u8], cigar: &mut Cigar) -> Result<(), ParseError> {
    let mut ops = Vec::from(mem::take(cigar));

    while !src.is_empty() {
        let op = parse_op(&mut src).map_err(ParseError::InvalidOp)?;
        ops.push(op);
    }

    *cigar = Cigar::try_from(ops).map_err(|_| ParseError::Invalid)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cigar() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::cigar::{op::Kind, Op};

        let src = b"1M13N144S";
        let mut cigar = Cigar::default();
        parse_cigar(src, &mut cigar)?;

        let expected = Cigar::try_from(vec![
            Op::new(Kind::Match, 1),
            Op::new(Kind::Skip, 13),
            Op::new(Kind::SoftClip, 144),
        ])?;

        assert_eq!(cigar, expected);

        Ok(())
    }
}
