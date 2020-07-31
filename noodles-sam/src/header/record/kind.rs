//! SAM header record kind.

use std::{error, fmt, str::FromStr};

/// A SAM header record kind.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kind {
    /// Header (`HD`).
    Header,
    /// Reference sequence (`SQ`).
    ReferenceSequence,
    /// Read group (`RG`).
    ReadGroup,
    /// Program (`PG`).
    Program,
    /// Comment (`CO`).
    Comment,
}

impl AsRef<str> for Kind {
    fn as_ref(&self) -> &str {
        match self {
            Self::Header => "HD",
            Self::ReferenceSequence => "SQ",
            Self::ReadGroup => "RG",
            Self::Program => "PG",
            Self::Comment => "CO",
        }
    }
}

impl fmt::Display for Kind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "@{}", self.as_ref())
    }
}

/// An error returned when a raw SAM header record kind fails to parse.
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

impl FromStr for Kind {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            "@HD" => Ok(Self::Header),
            "@SQ" => Ok(Self::ReferenceSequence),
            "@RG" => Ok(Self::ReadGroup),
            "@PG" => Ok(Self::Program),
            "@CO" => Ok(Self::Comment),
            _ => Err(ParseError::Invalid),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert_eq!("@HD".parse(), Ok(Kind::Header));
        assert_eq!("@SQ".parse(), Ok(Kind::ReferenceSequence));
        assert_eq!("@RG".parse(), Ok(Kind::ReadGroup));
        assert_eq!("@PG".parse(), Ok(Kind::Program));
        assert_eq!("@CO".parse(), Ok(Kind::Comment));

        assert_eq!("".parse::<Kind>(), Err(ParseError::Empty));
        assert_eq!("@NO".parse::<Kind>(), Err(ParseError::Invalid));
        assert_eq!("HD".parse::<Kind>(), Err(ParseError::Invalid));
    }
}
