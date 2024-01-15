//! SAM header reference sequence name.

use std::{borrow::Borrow, error, fmt, ops::Deref, str::FromStr};

/// A SAM header reference sequence name.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Name(String);

impl Borrow<str> for Name {
    fn borrow(&self) -> &str {
        &self.0
    }
}

impl Borrow<[u8]> for Name {
    fn borrow(&self) -> &[u8] {
        self.0.as_bytes()
    }
}

impl Deref for Name {
    type Target = String;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for Name {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// An error returned when a raw SAM record reference sequence name fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid(String),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid(s) => write!(f, "invalid input: {s}"),
        }
    }
}

impl FromStr for Name {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            Err(ParseError::Empty)
        } else if !is_valid_name(s) {
            Err(ParseError::Invalid(s.into()))
        } else {
            Ok(Self(s.into()))
        }
    }
}

// § 1.2.1 Character set restrictions (2021-01-07)
fn is_valid_name_char(c: char) -> bool {
    ('!'..='~').contains(&c)
        && !matches!(
            c,
            '\\' | ',' | '"' | '`' | '\'' | '(' | ')' | '[' | ']' | '{' | '}' | '<' | '>',
        )
}

pub(crate) fn is_valid_name(s: &str) -> bool {
    let mut chars = s.chars();

    if let Some(c) = chars.next() {
        if c == '*' || c == '=' || !is_valid_name_char(c) {
            return false;
        }

        chars.all(is_valid_name_char)
    } else {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() -> Result<(), ParseError> {
        let reference_sequence_name: Name = "sq0".parse()?;
        assert_eq!(reference_sequence_name.to_string(), "sq0");

        Ok(())
    }

    #[test]
    fn test_from_str() {
        assert_eq!("sq0".parse(), Ok(Name(String::from("sq0"))));
        assert_eq!("sq0*".parse(), Ok(Name(String::from("sq0*"))));
        assert_eq!("sq0=".parse(), Ok(Name(String::from("sq0="))));

        assert_eq!("".parse::<Name>(), Err(ParseError::Empty));

        assert_eq!(
            "sq 0".parse::<Name>(),
            Err(ParseError::Invalid(String::from("sq 0")))
        );

        assert_eq!(
            "sq[0]".parse::<Name>(),
            Err(ParseError::Invalid(String::from("sq[0]")))
        );

        assert_eq!(
            ">sq0".parse::<Name>(),
            Err(ParseError::Invalid(String::from(">sq0")))
        );

        assert_eq!(
            "*sq0".parse::<Name>(),
            Err(ParseError::Invalid(String::from("*sq0")))
        );

        assert_eq!(
            "=sq0".parse::<Name>(),
            Err(ParseError::Invalid(String::from("=sq0")))
        );
    }
}
