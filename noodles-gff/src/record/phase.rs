//! GFF record phase.

use std::{error, fmt, str::FromStr};

/// A GFF record phase.
///
/// The phase is used for CDS (coding sequence) features to describe where the next codon begins
/// relative to the 5' end.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Phase {
    /// The codon begins at the first nucleotide (`0`).
    Zero,
    /// The codon begins at the second nucleotide (`1`).
    One,
    /// The codon begins at the third nucleotide (`2`).
    Two,
}

impl AsRef<str> for Phase {
    fn as_ref(&self) -> &str {
        match self {
            Self::Zero => "0",
            Self::One => "1",
            Self::Two => "2",
        }
    }
}

impl fmt::Display for Phase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

/// An error returned when a raw GFF record phase fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The phase is invalid.
    Invalid(String),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid(s) => write!(f, "expected {{0, 1, 2}}, got {s}"),
        }
    }
}

impl FromStr for Phase {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            "0" => Ok(Self::Zero),
            "1" => Ok(Self::One),
            "2" => Ok(Self::Two),
            _ => Err(ParseError::Invalid(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Phase::Zero.to_string(), "0");
        assert_eq!(Phase::One.to_string(), "1");
        assert_eq!(Phase::Two.to_string(), "2");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("0".parse::<Phase>()?, Phase::Zero);
        assert_eq!("1".parse::<Phase>()?, Phase::One);
        assert_eq!("2".parse::<Phase>()?, Phase::Two);

        assert_eq!("".parse::<Phase>(), Err(ParseError::Empty));
        assert_eq!(
            "3".parse::<Phase>(),
            Err(ParseError::Invalid(String::from("3")))
        );

        Ok(())
    }
}
