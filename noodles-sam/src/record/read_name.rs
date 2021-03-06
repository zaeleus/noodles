//! SAM record read name.

use std::{error, fmt, ops::Deref, str::FromStr};

// § 1.4 The alignment section: mandatory fields (2020-07-19)
const MAX_LENGTH: usize = 254;

/// A SAM record read name.
///
/// This is also called a query name.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReadName(String);

impl Deref for ReadName {
    type Target = String;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for ReadName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

/// An error returned when a raw SAM record read name fails to parse.
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

impl FromStr for ReadName {
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

// § 1.2.1 Character set restrictions (2021-01-07)
fn is_valid_name_char(c: char) -> bool {
    ('!'..='~').contains(&c) && c != '@'
}

fn is_valid_name(s: &str) -> bool {
    s.len() <= MAX_LENGTH && s.chars().all(is_valid_name_char)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() -> Result<(), ParseError> {
        let read_name: ReadName = "r0".parse()?;
        assert_eq!(read_name.to_string(), "r0");
        Ok(())
    }

    #[test]
    fn test_from_str() {
        assert_eq!("r0".parse(), Ok(ReadName(String::from("r0"))));

        assert_eq!("".parse::<ReadName>(), Err(ParseError::Empty));
        assert_eq!("r 0".parse::<ReadName>(), Err(ParseError::Invalid));
        assert_eq!("@r0".parse::<ReadName>(), Err(ParseError::Invalid));

        let s: String = (0..MAX_LENGTH + 1).map(|_| 'N').collect();
        assert_eq!(s.parse::<ReadName>(), Err(ParseError::Invalid));
    }
}
