//! VCF header header contig record key.

use std::{error, fmt, str::FromStr};

/// A VCF header header contig record key.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Key {
    /// (`ID`).
    Id,
    /// (`length`).
    Length,
    /// (`IDX`).
    Idx,
    /// Any other key.
    Other(String),
}

impl AsRef<str> for Key {
    fn as_ref(&self) -> &str {
        match self {
            Self::Id => "ID",
            Self::Length => "length",
            Self::Idx => "IDX",
            Self::Other(s) => s,
        }
    }
}

impl fmt::Display for Key {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

/// An error returned when a raw VCF header contig record key fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
        }
    }
}

impl FromStr for Key {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            "ID" => Ok(Self::Id),
            "length" => Ok(Self::Length),
            "IDX" => Ok(Self::Idx),
            _ => Ok(Self::Other(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Key::Id.to_string(), "ID");
        assert_eq!(Key::Length.to_string(), "length");
        assert_eq!(Key::Idx.to_string(), "IDX");
        assert_eq!(Key::Other(String::from("md5")).to_string(), "md5");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("ID".parse(), Ok(Key::Id));
        assert_eq!("length".parse(), Ok(Key::Length));
        assert_eq!("IDX".parse(), Ok(Key::Idx));
        assert_eq!("assembly".parse(), Ok(Key::Other(String::from("assembly"))));
        assert_eq!("md5".parse(), Ok(Key::Other(String::from("md5"))));
        assert_eq!("species".parse(), Ok(Key::Other(String::from("species"))));
        assert_eq!("taxonomy".parse(), Ok(Key::Other(String::from("taxonomy"))));
        assert_eq!("URL".parse(), Ok(Key::Other(String::from("URL"))));

        assert_eq!("".parse::<Key>(), Err(ParseError::Empty));
    }
}
