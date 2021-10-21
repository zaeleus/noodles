//! SAM record sequence and bases.

mod base;

pub use self::base::Base;

use std::{error, fmt, ops::Deref, str::FromStr};

use super::NULL_FIELD;

/// A SAM record sequence.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Sequence(Vec<Base>);

impl Deref for Sequence {
    type Target = [Base];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<Vec<Base>> for Sequence {
    fn from(bases: Vec<Base>) -> Self {
        Self(bases)
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.is_empty() {
            write!(f, "{}", NULL_FIELD)
        } else {
            for base in &self.0 {
                write!(f, "{}", base)?;
            }

            Ok(())
        }
    }
}

/// An error returned when a raw SAM record sequence fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The raw sequence has an invalid base.
    InvalidBase(base::TryFromCharError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::InvalidBase(e) => write!(f, "invalid base: {}", e),
        }
    }
}

impl FromStr for Sequence {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            NULL_FIELD => Ok(Self::default()),
            _ => s
                .chars()
                .map(|c| c.to_ascii_uppercase())
                .map(Base::try_from)
                .collect::<Result<Vec<_>, _>>()
                .map(Self::from)
                .map_err(ParseError::InvalidBase),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let sequence = Sequence::from(vec![Base::A, Base::T, Base::C, Base::G]);
        assert_eq!(sequence.to_string(), "ATCG");
    }

    #[test]
    fn test_from_str() {
        let expected = Sequence(vec![Base::A, Base::T, Base::C, Base::G]);
        assert_eq!("ATCG".parse::<Sequence>(), Ok(expected));

        let expected = Sequence(vec![Base::A, Base::T, Base::C, Base::G]);
        assert_eq!("atcg".parse::<Sequence>(), Ok(expected));

        let expected = Sequence(vec![Base::A, Base::T, Base::C, Base::G]);
        assert_eq!("aTcG".parse::<Sequence>(), Ok(expected));

        assert_eq!("*".parse::<Sequence>(), Ok(Sequence::default()));

        assert_eq!("".parse::<Sequence>(), Err(ParseError::Empty));
    }
}
