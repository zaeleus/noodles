//! SAM header program tag.

use std::{error, fmt, str::FromStr};

/// A SAM header program tag.
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub enum Tag {
    /// Program ID (`ID`).
    Id,
    /// Program name (`PN`).
    Name,
    /// Command line (`CL`).
    CommandLine,
    /// Previous program ID (`PP`).
    PreviousId,
    /// Description (`DS`).
    Description,
    /// Program version (`VN`).
    Version,
    /// Any other program tag.
    Other(String),
}

impl AsRef<str> for Tag {
    fn as_ref(&self) -> &str {
        match self {
            Self::Id => "ID",
            Self::Name => "PN",
            Self::CommandLine => "CL",
            Self::PreviousId => "PP",
            Self::Description => "DS",
            Self::Version => "VN",
            Self::Other(s) => s,
        }
    }
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_ref())
    }
}

/// An error returned when a raw SAM header program tag fails to parse.
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

impl FromStr for Tag {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            "ID" => Ok(Self::Id),
            "PN" => Ok(Self::Name),
            "CL" => Ok(Self::CommandLine),
            "PP" => Ok(Self::PreviousId),
            "DS" => Ok(Self::Description),
            "VN" => Ok(Self::Version),
            _ => {
                if s.len() == 2 {
                    Ok(Self::Other(s.into()))
                } else {
                    Err(ParseError::Invalid)
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Tag::Id.to_string(), "ID");
        assert_eq!(Tag::Name.to_string(), "PN");
        assert_eq!(Tag::CommandLine.to_string(), "CL");
        assert_eq!(Tag::PreviousId.to_string(), "PP");
        assert_eq!(Tag::Description.to_string(), "DS");
        assert_eq!(Tag::Version.to_string(), "VN");
        assert_eq!(Tag::Other(String::from("ND")).to_string(), "ND");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("ID".parse(), Ok(Tag::Id));
        assert_eq!("PN".parse(), Ok(Tag::Name));
        assert_eq!("CL".parse(), Ok(Tag::CommandLine));
        assert_eq!("PP".parse(), Ok(Tag::PreviousId));
        assert_eq!("DS".parse(), Ok(Tag::Description));
        assert_eq!("VN".parse(), Ok(Tag::Version));
        assert_eq!("ND".parse(), Ok(Tag::Other(String::from("ND"))));

        assert_eq!("".parse::<Tag>(), Err(ParseError::Empty));
        assert_eq!("NDL".parse::<Tag>(), Err(ParseError::Invalid));
    }
}
