pub(crate) mod op;

use std::{error, fmt};

use self::op::parse_op;
use crate::alignment::record_buf::Cigar;

/// An error when a raw SAM record CIGAR fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
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
            Self::InvalidOp(_) => write!(f, "invalid op"),
        }
    }
}

pub(super) fn parse_cigar(mut src: &[u8], cigar: &mut Cigar) -> Result<(), ParseError> {
    if src.is_empty() {
        return Err(ParseError::Empty);
    }

    cigar.as_mut().clear();

    while !src.is_empty() {
        let op = parse_op(&mut src).map_err(ParseError::InvalidOp)?;
        cigar.as_mut().push(op);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cigar() -> Result<(), ParseError> {
        use crate::alignment::record::cigar::{op::Kind, Op};

        let mut cigar = Cigar::default();

        cigar.as_mut().clear();
        parse_cigar(b"8M13N", &mut cigar)?;
        assert_eq!(
            cigar,
            [Op::new(Kind::Match, 8), Op::new(Kind::Skip, 13),]
                .into_iter()
                .collect()
        );

        cigar.as_mut().clear();
        parse_cigar(b"8M13N144S", &mut cigar)?;
        assert_eq!(
            cigar,
            [
                Op::new(Kind::Match, 8),
                Op::new(Kind::Skip, 13),
                Op::new(Kind::SoftClip, 144),
            ]
            .into_iter()
            .collect()
        );

        cigar.as_mut().clear();
        assert_eq!(parse_cigar(b"", &mut cigar), Err(ParseError::Empty));

        cigar.as_mut().clear();
        assert!(matches!(
            parse_cigar(b"8Z", &mut cigar),
            Err(ParseError::InvalidOp(_))
        ));

        Ok(())
    }
}
