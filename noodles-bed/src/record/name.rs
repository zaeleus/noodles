//! BED record name.

use std::{error, fmt, ops::Deref, str::FromStr};

const MAX_LENGTH: usize = 255;

/// A BED record name.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Name(String);

impl Deref for Name {
    type Target = str;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for Name {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self)
    }
}

/// An error returned when a raw BED record name fails to parse.
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

impl FromStr for Name {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            Err(ParseError::Empty)
        } else if !is_valid_name(s) {
            Err(ParseError::Invalid)
        } else {
            Ok(Self(s.into()))
        }
    }
}

fn is_valid_name_char(c: char) -> bool {
    matches!(c, ' '..='~')
}

fn is_valid_name(s: &str) -> bool {
    s.len() <= MAX_LENGTH && s.chars().all(is_valid_name_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let name = Name(String::from("ndls1"));
        assert_eq!(name.to_string(), "ndls1");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("ndls1".parse(), Ok(Name(String::from("ndls1"))));
        assert_eq!(" ~".parse(), Ok(Name(String::from(" ~"))));

        assert_eq!("".parse::<Name>(), Err(ParseError::Empty));
        assert_eq!("🍜".parse::<Name>(), Err(ParseError::Invalid));

        let s = "n".repeat(MAX_LENGTH + 1);
        assert_eq!(s.parse::<Name>(), Err(ParseError::Invalid));
    }
}
