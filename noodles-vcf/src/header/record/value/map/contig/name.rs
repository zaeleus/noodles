//! VCF header contig name.

use std::{error, fmt, str::FromStr};

/// A VCF header contig name.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Name(String);

impl AsRef<str> for Name {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for Name {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// An error returned when a raw VCF header contig name fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => "invalid input".fmt(f),
        }
    }
}

impl FromStr for Name {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use crate::record::chromosome::is_valid_name;

        if is_valid_name(s) {
            Ok(Self(s.into()))
        } else {
            Err(ParseError::Invalid)
        }
    }
}
