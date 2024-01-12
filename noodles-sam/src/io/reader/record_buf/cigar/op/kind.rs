use std::{error, fmt};

use crate::alignment::record::cigar::op::Kind;

/// An error when a raw SAM record CIGAR op kind fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// The input is invalid.
    Invalid { actual: u8 },
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::Invalid { actual } => write!(
                f,
                "invalid input: expected {{M, I, D, N, S, H, P, =, X}}, got {c}",
                c = char::from(*actual)
            ),
        }
    }
}

pub(super) fn parse_kind(src: &mut &[u8]) -> Result<Kind, ParseError> {
    let (n, rest) = src.split_first().ok_or(ParseError::UnexpectedEof)?;

    *src = rest;

    match n {
        b'M' => Ok(Kind::Match),
        b'I' => Ok(Kind::Insertion),
        b'D' => Ok(Kind::Deletion),
        b'N' => Ok(Kind::Skip),
        b'S' => Ok(Kind::SoftClip),
        b'H' => Ok(Kind::HardClip),
        b'P' => Ok(Kind::Pad),
        b'=' => Ok(Kind::SequenceMatch),
        b'X' => Ok(Kind::SequenceMismatch),
        _ => Err(ParseError::Invalid { actual: *n }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_kind() {
        fn t(mut src: &[u8], expected: Kind) {
            assert_eq!(parse_kind(&mut src), Ok(expected));
        }

        t(b"M", Kind::Match);
        t(b"I", Kind::Insertion);
        t(b"D", Kind::Deletion);
        t(b"N", Kind::Skip);
        t(b"S", Kind::SoftClip);
        t(b"H", Kind::HardClip);
        t(b"P", Kind::Pad);
        t(b"=", Kind::SequenceMatch);
        t(b"X", Kind::SequenceMismatch);

        assert_eq!(parse_kind(&mut &[][..]), Err(ParseError::UnexpectedEof));

        assert_eq!(
            parse_kind(&mut &b"!"[..]),
            Err(ParseError::Invalid { actual: b'!' })
        );
    }
}
