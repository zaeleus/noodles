//! SAM header program tag.

use std::{
    error,
    fmt::{self, Write},
    str::FromStr,
};

const LENGTH: usize = 2;

/// A SAM header program tag.
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
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
    Other([u8; LENGTH]),
}

impl AsRef<[u8; LENGTH]> for Tag {
    fn as_ref(&self) -> &[u8; LENGTH] {
        match self {
            Self::Id => b"ID",
            Self::Name => b"PN",
            Self::CommandLine => b"CL",
            Self::PreviousId => b"PP",
            Self::Description => b"DS",
            Self::Version => b"VN",
            Self::Other(tag) => tag,
        }
    }
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let bytes = self.as_ref();
        f.write_char(char::from(bytes[0]))?;
        f.write_char(char::from(bytes[1]))?;
        Ok(())
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
        let bytes = s.as_bytes();

        match bytes.len() {
            0 => Err(ParseError::Empty),
            2 => Self::try_from([bytes[0], bytes[1]]),
            _ => Err(ParseError::Invalid),
        }
    }
}

impl TryFrom<[u8; LENGTH]> for Tag {
    type Error = ParseError;

    fn try_from(b: [u8; LENGTH]) -> Result<Self, Self::Error> {
        match &b {
            b"ID" => Ok(Self::Id),
            b"PN" => Ok(Self::Name),
            b"CL" => Ok(Self::CommandLine),
            b"PP" => Ok(Self::PreviousId),
            b"DS" => Ok(Self::Description),
            b"VN" => Ok(Self::Version),
            _ => Ok(Self::Other(b)),
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
        assert_eq!(Tag::Other([b'N', b'D']).to_string(), "ND");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("ID".parse(), Ok(Tag::Id));
        assert_eq!("PN".parse(), Ok(Tag::Name));
        assert_eq!("CL".parse(), Ok(Tag::CommandLine));
        assert_eq!("PP".parse(), Ok(Tag::PreviousId));
        assert_eq!("DS".parse(), Ok(Tag::Description));
        assert_eq!("VN".parse(), Ok(Tag::Version));
        assert_eq!("ND".parse(), Ok(Tag::Other([b'N', b'D'])));

        assert_eq!("".parse::<Tag>(), Err(ParseError::Empty));
        assert_eq!("NDL".parse::<Tag>(), Err(ParseError::Invalid));
    }
}
