//! VCF header record other key.

use std::{borrow::Borrow, error, fmt, str::FromStr};

/// A nonstandard VCF record key.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Other(pub(super) String);

impl AsRef<str> for Other {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

impl Borrow<str> for Other {
    fn borrow(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for Other {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.as_ref().fmt(f)
    }
}

/// An error returned when a raw VCF header record other key fails to a parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => write!(f, "invalid input"),
        }
    }
}

impl FromStr for Other {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use super::Key;

        match Key::from(s) {
            Key::Standard(_) => Err(ParseError::Invalid),
            Key::Other(key) => Ok(key),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert_eq!("NOODLES".parse(), Ok(Other(String::from("NOODLES"))));
        assert_eq!("INFO".parse::<Other>(), Err(ParseError::Invalid));
    }
}
