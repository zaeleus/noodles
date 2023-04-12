use std::{error, fmt, mem};

use crate::record::{
    sequence::{base, Base},
    Sequence,
};

/// An error when a raw SAM record sequence fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A base is invalid.
    InvalidBase(base::TryFromCharError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidBase(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBase(_) => write!(f, "invalid base"),
        }
    }
}

pub(crate) fn parse_sequence(src: &[u8], sequence: &mut Sequence) -> Result<(), ParseError> {
    let mut bases = Vec::from(mem::take(sequence));

    for &n in src {
        let base = Base::try_from(n).map_err(ParseError::InvalidBase)?;
        bases.push(base);
    }

    *sequence = Sequence::from(bases);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_sequence() -> Result<(), ParseError> {
        let mut sequence = Sequence::default();

        sequence.clear();
        parse_sequence(b"", &mut sequence)?;
        let expected = Sequence::default();
        assert_eq!(sequence, expected);

        sequence.clear();
        parse_sequence(b"ACGT", &mut sequence)?;
        let expected = Sequence::from(vec![Base::A, Base::C, Base::G, Base::T]);
        assert_eq!(sequence, expected);

        sequence.clear();
        assert!(matches!(
            parse_sequence(&[0x07], &mut sequence),
            Err(ParseError::InvalidBase(_))
        ));

        Ok(())
    }
}
