//! VCF record genotype value allele phasing.

use std::{error, fmt, str::FromStr};

/// A VCF record genotype value allele phasing.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Phasing {
    /// The allele is phased.
    Phased,
    /// The allele is unphased.
    Unphased,
}

impl AsRef<str> for Phasing {
    fn as_ref(&self) -> &str {
        match self {
            Self::Phased => "|",
            Self::Unphased => "/",
        }
    }
}

impl fmt::Display for Phasing {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

/// An error returned when a raw VCF record genotype value allele phasing fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid => f.write_str("invalid input"),
        }
    }
}

impl FromStr for Phasing {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            "|" => Ok(Self::Phased),
            "/" => Ok(Self::Unphased),
            _ => Err(ParseError::Invalid),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Phasing::Phased.to_string(), "|");
        assert_eq!(Phasing::Unphased.to_string(), "/");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("|".parse(), Ok(Phasing::Phased));
        assert_eq!("/".parse(), Ok(Phasing::Unphased));

        assert_eq!("".parse::<Phasing>(), Err(ParseError::Empty));
        assert_eq!(":".parse::<Phasing>(), Err(ParseError::Invalid));
    }
}
