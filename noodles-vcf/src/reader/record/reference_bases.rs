mod base;

use std::{error, fmt};

use self::base::parse_base;
use crate::record::ReferenceBases;

/// An error when raw VCF record reference bases fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A base is invalid.
    InvalidBase(base::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidBase(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::InvalidBase(_) => write!(f, "invalid base"),
        }
    }
}

pub(super) fn parse_reference_bases(
    s: &str,
    reference_bases: &mut ReferenceBases,
) -> Result<(), ParseError> {
    if s.is_empty() {
        return Err(ParseError::Empty);
    }

    reference_bases.0.clear();

    for c in s.chars() {
        let base = parse_base(c).map_err(ParseError::InvalidBase)?;
        reference_bases.0.push(base);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_reference_bases() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::reference_bases::Base;

        let mut reference_bases = ReferenceBases::try_from(vec![Base::N])?;

        let expected = [Base::A, Base::T, Base::C, Base::G, Base::N];

        parse_reference_bases("ATCGN", &mut reference_bases)?;
        assert_eq!(&reference_bases[..], &expected[..]);

        parse_reference_bases("atcgn", &mut reference_bases)?;
        assert_eq!(&reference_bases[..], &expected[..]);

        parse_reference_bases("AtCgN", &mut reference_bases)?;
        assert_eq!(&reference_bases[..], &expected[..]);

        assert_eq!(
            parse_reference_bases("", &mut reference_bases),
            Err(ParseError::Empty)
        );

        assert!(matches!(
            parse_reference_bases(".", &mut reference_bases),
            Err(ParseError::InvalidBase(_))
        ));

        assert!(matches!(
            parse_reference_bases("Z", &mut reference_bases),
            Err(ParseError::InvalidBase(_))
        ));

        Ok(())
    }
}
