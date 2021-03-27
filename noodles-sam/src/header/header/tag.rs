//! SAM header header tag.

use std::{error, fmt, str::FromStr};

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
    Other(String),
}

impl AsRef<str> for Tag {
    fn as_ref(&self) -> &str {
        match self {
            Self::Version => "VN",
            Self::SortOrder => "SO",
            Self::GroupOrder => "GO",
            Self::SubsortOrder => "SS",
            Self::Other(s) => s,
        }
    }
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_ref())
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
        match s {
            "" => Err(ParseError::Empty),
            "VN" => Ok(Self::Version),
            "SO" => Ok(Self::SortOrder),
            "GO" => Ok(Self::GroupOrder),
            "SS" => Ok(Self::SubsortOrder),
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
        assert_eq!(Tag::Version.to_string(), "VN");
        assert_eq!(Tag::SortOrder.to_string(), "SO");
        assert_eq!(Tag::GroupOrder.to_string(), "GO");
        assert_eq!(Tag::SubsortOrder.to_string(), "SS");
        assert_eq!(Tag::Other(String::from("ND")).to_string(), "ND");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("VN".parse::<Tag>(), Ok(Tag::Version));
        assert_eq!("SO".parse::<Tag>(), Ok(Tag::SortOrder));
        assert_eq!("GO".parse::<Tag>(), Ok(Tag::GroupOrder));
        assert_eq!("SS".parse::<Tag>(), Ok(Tag::SubsortOrder));
        assert_eq!("ND".parse::<Tag>(), Ok(Tag::Other(String::from("ND"))));

        assert_eq!("".parse::<Tag>(), Err(ParseError::Empty));
        assert_eq!("NDL".parse::<Tag>(), Err(ParseError::Invalid));
    }
}
