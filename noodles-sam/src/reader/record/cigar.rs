mod op;

use std::{error, fmt, mem};

use self::op::parse_op;
use crate::record::Cigar;

/// An error when a raw SAM record CIGAR fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
    /// An op is invalid.
    InvalidOp(op::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidOp(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::Invalid => write!(f, "invalid input"),
            Self::InvalidOp(_) => write!(f, "invalid op"),
        }
    }
}

pub(crate) fn parse_cigar(mut src: &[u8], cigar: &mut Cigar) -> Result<(), ParseError> {
    if src.is_empty() {
        return Err(ParseError::Empty);
    }

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

        let mut cigar = Cigar::default();

        cigar.clear();
        parse_cigar(b"8M13N", &mut cigar)?;
        assert_eq!(
            cigar,
            Cigar::try_from(vec![Op::new(Kind::Match, 8), Op::new(Kind::Skip, 13),])?
        );

        cigar.clear();
        parse_cigar(b"8M13N144S", &mut cigar)?;
        assert_eq!(
            cigar,
            Cigar::try_from(vec![
                Op::new(Kind::Match, 8),
                Op::new(Kind::Skip, 13),
                Op::new(Kind::SoftClip, 144),
            ])?
        );

        cigar.clear();
        assert_eq!(parse_cigar(b"", &mut cigar), Err(ParseError::Empty));

        cigar.clear();
        assert!(matches!(
            parse_cigar(b"8Z", &mut cigar),
            Err(ParseError::InvalidOp(_))
        ));

        Ok(())
    }
}
