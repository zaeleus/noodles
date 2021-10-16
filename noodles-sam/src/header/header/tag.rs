//! SAM header header tag.

use std::{
    convert::TryFrom,
    error,
    fmt::{self, Write},
    str::FromStr,
};

const LENGTH: usize = 2;

/// A SAM header header tag.
#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub enum Tag {
    /// Format version (`VN`).
    Version,
    /// Sorting order of alignments (`SO`).
    SortOrder,
    /// Group order of alignments (`GO`).
    GroupOrder,
    /// Subsort order of alignments (`SS`).
    SubsortOrder,
    /// Any other header tag.
    Other([u8; LENGTH]),
}

impl AsRef<[u8; LENGTH]> for Tag {
    fn as_ref(&self) -> &[u8; LENGTH] {
        match self {
            Self::Version => b"VN",
            Self::SortOrder => b"SO",
            Self::GroupOrder => b"GO",
            Self::SubsortOrder => b"SS",
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

/// An error returned when a raw SAM header header tag fails to parse.
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
            b"VN" => Ok(Self::Version),
            b"SO" => Ok(Self::SortOrder),
            b"GO" => Ok(Self::GroupOrder),
            b"SS" => Ok(Self::SubsortOrder),
            _ => Ok(Self::Other(b)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Tag::Version.to_string(), "VN");
        assert_eq!(Tag::SortOrder.to_string(), "SO");
        assert_eq!(Tag::GroupOrder.to_string(), "GO");
        assert_eq!(Tag::SubsortOrder.to_string(), "SS");
        assert_eq!(Tag::Other([b'N', b'D']).to_string(), "ND");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("VN".parse(), Ok(Tag::Version));
        assert_eq!("SO".parse(), Ok(Tag::SortOrder));
        assert_eq!("GO".parse(), Ok(Tag::GroupOrder));
        assert_eq!("SS".parse(), Ok(Tag::SubsortOrder));
        assert_eq!("ND".parse(), Ok(Tag::Other([b'N', b'D'])));

        assert_eq!("".parse::<Tag>(), Err(ParseError::Empty));
        assert_eq!("NDL".parse::<Tag>(), Err(ParseError::Invalid));
    }
}
