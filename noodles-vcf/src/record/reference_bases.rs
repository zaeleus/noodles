mod base;

pub use self::base::Base;

use std::{convert::TryFrom, error, fmt, ops::Deref, str::FromStr};

use super::MISSING_FIELD;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReferenceBases(Vec<Base>);

impl Deref for ReferenceBases {
    type Target = [Base];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for ReferenceBases {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for base in self.iter() {
            write!(f, "{}", char::from(*base))?;
        }

        Ok(())
    }
}

/// An error returned when raw VCF reference bases fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is missing (`.`).
    Missing,
    /// The input has an invalid base.
    InvalidBase(base::TryFromCharError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Missing => f.write_str("missing input"),
            Self::InvalidBase(e) => write!(f, "invalid base: {}", e),
        }
    }
}

impl FromStr for ReferenceBases {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            MISSING_FIELD => Err(ParseError::Missing),
            _ => s
                .chars()
                .map(|c| c.to_ascii_uppercase())
                .map(Base::try_from)
                .collect::<Result<_, _>>()
                .map(ReferenceBases)
                .map_err(ParseError::InvalidBase),
        }
    }
}

/// An error returned when a vector of bases fails to convert to reference bases.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromBaseVectorError {
    /// The input is empty.
    Empty,
}

impl error::Error for TryFromBaseVectorError {}

impl fmt::Display for TryFromBaseVectorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty list"),
        }
    }
}

impl TryFrom<Vec<Base>> for ReferenceBases {
    type Error = TryFromBaseVectorError;

    fn try_from(bases: Vec<Base>) -> Result<Self, Self::Error> {
        if bases.is_empty() {
            Err(TryFromBaseVectorError::Empty)
        } else {
            Ok(Self(bases))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let reference_bases = ReferenceBases(vec![Base::A, Base::T, Base::C, Base::G, Base::N]);
        assert_eq!(reference_bases.to_string(), "ATCGN");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let expected = [Base::A, Base::T, Base::C, Base::G, Base::N];

        let bases: ReferenceBases = "ATCGN".parse()?;
        assert_eq!(&bases[..], &expected[..]);

        let bases: ReferenceBases = "atcgn".parse()?;
        assert_eq!(&bases[..], &expected[..]);

        let bases: ReferenceBases = "AtCgN".parse()?;
        assert_eq!(&bases[..], &expected[..]);

        assert_eq!("".parse::<ReferenceBases>(), Err(ParseError::Empty));
        assert_eq!(".".parse::<ReferenceBases>(), Err(ParseError::Missing));
        assert!(matches!(
            "Z".parse::<ReferenceBases>(),
            Err(ParseError::InvalidBase(_))
        ));

        Ok(())
    }

    #[test]
    fn test_try_from_vec_bases_for_reference_bases() {
        assert_eq!(
            ReferenceBases::try_from(vec![Base::A]),
            Ok(ReferenceBases(vec![Base::A]))
        );
        assert_eq!(
            ReferenceBases::try_from(Vec::new()),
            Err(TryFromBaseVectorError::Empty)
        );
    }
}
