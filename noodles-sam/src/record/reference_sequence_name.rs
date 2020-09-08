//! SAM record reference sequence name.

use std::{error, fmt, ops::Deref, str::FromStr};

/// A SAM record reference sequence name.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceSequenceName(String);

impl Deref for ReferenceSequenceName {
    type Target = String;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for ReferenceSequenceName {
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

impl FromStr for ReferenceSequenceName {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            Err(ParseError::Empty)
        } else if s.starts_with('*') || s.starts_with('=') {
            Err(ParseError::Invalid)
        } else {
            Ok(ReferenceSequenceName(s.into()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() -> Result<(), ParseError> {
        let reference_sequence_name: ReferenceSequenceName = "sq0".parse()?;
        assert_eq!(reference_sequence_name.to_string(), "sq0");

        Ok(())
    }

    #[test]
    fn test_from_str() {
        assert_eq!(
            "sq0".parse(),
            Ok(ReferenceSequenceName(String::from("sq0")))
        );

        assert_eq!(
            "sq0*".parse(),
            Ok(ReferenceSequenceName(String::from("sq0*")))
        );

        assert_eq!(
            "sq0=".parse(),
            Ok(ReferenceSequenceName(String::from("sq0=")))
        );

        assert_eq!("".parse::<ReferenceSequenceName>(), Err(ParseError::Empty));

        assert_eq!(
            "*sq0".parse::<ReferenceSequenceName>(),
            Err(ParseError::Invalid)
        );

        assert_eq!(
            "=sq0".parse::<ReferenceSequenceName>(),
            Err(ParseError::Invalid)
        );

        assert_eq!(
            "*".parse::<ReferenceSequenceName>(),
            Err(ParseError::Invalid)
        );

        assert_eq!(
            "=".parse::<ReferenceSequenceName>(),
            Err(ParseError::Invalid)
        );
    }
}
