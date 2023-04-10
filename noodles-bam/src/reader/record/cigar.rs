mod op;

use std::{error, fmt, mem};

use bytes::Buf;
use noodles_sam::record::Cigar;

pub(super) use self::op::decode_op;

/// An error when a raw BAM record CIGAR fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// An op is invalid.
    InvalidOp(op::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::UnexpectedEof => None,
            Self::InvalidOp(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::InvalidOp(_) => write!(f, "invalid op"),
        }
    }
}

pub fn get_cigar<B>(src: &mut B, cigar: &mut Cigar, n_cigar_op: usize) -> Result<(), ParseError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u32>() * n_cigar_op {
        return Err(ParseError::UnexpectedEof);
    }

    cigar.clear();

    for _ in 0..n_cigar_op {
        let op = decode_op(src.get_u32_le()).map_err(ParseError::InvalidOp)?;
        cigar.as_mut().push(op);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_cigar() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            mut src: &[u8],
            actual: &mut Cigar,
            n_cigar_op: usize,
            expected: &Cigar,
        ) -> Result<(), ParseError> {
            get_cigar(&mut src, actual, n_cigar_op)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        let mut buf = Cigar::default();

        t(&[], &mut buf, 0, &Cigar::default())?;
        t(&[0x40, 0x00, 0x00, 0x00], &mut buf, 1, &"4M".parse()?)?;
        t(
            &[0x40, 0x00, 0x00, 0x00, 0x25, 0x00, 0x00, 0x00],
            &mut buf,
            2,
            &"4M2H".parse()?,
        )?;

        Ok(())
    }
}
