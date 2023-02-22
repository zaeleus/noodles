//! VCF record ID.

use std::{error, fmt, ops::Deref, str::FromStr};

/// A VCF record ID.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Id(String);

impl Deref for Id {
    type Target = str;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for Id {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

/// An error returned when a raw VCF record ID fails to parse.
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

impl FromStr for Id {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            Err(ParseError::Empty)
        } else if is_valid_id(s) {
            Ok(Self(s.into()))
        } else {
            Err(ParseError::Invalid)
        }
    }
}

fn is_valid_id(s: &str) -> bool {
    const MISSING: &str = ".";

    !s.is_empty() && s != MISSING && s.chars().all(|c| !c.is_ascii_whitespace())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert_eq!("nd0".parse(), Ok(Id(String::from("nd0"))));

        assert_eq!("".parse::<Id>(), Err(ParseError::Empty));
        assert_eq!(".".parse::<Id>(), Err(ParseError::Invalid));
        assert_eq!("nd 0".parse::<Id>(), Err(ParseError::Invalid));
    }
}
