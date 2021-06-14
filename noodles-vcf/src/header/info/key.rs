//! VCF header information record key.

use std::{error, fmt, str::FromStr};

/// A VCF header information record key.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Key {
    /// (`ID`).
    Id,
    /// (`Number`).
    Number,
    /// (`Type`).
    Type,
    /// (`Description`).
    Description,
    /// (`IDX`).
    Idx,
}

impl AsRef<str> for Key {
    fn as_ref(&self) -> &str {
        match self {
            Self::Id => "ID",
            Self::Number => "Number",
            Self::Type => "Type",
            Self::Description => "Description",
            Self::Idx => "IDX",
        }
    }
}

impl fmt::Display for Key {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

/// An error returned when a raw VCF header information record key fails to parse.
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

impl FromStr for Key {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            "ID" => Ok(Self::Id),
            "Number" => Ok(Self::Number),
            "Type" => Ok(Self::Type),
            "Description" => Ok(Self::Description),
            "IDX" => Ok(Self::Idx),
            _ => Err(ParseError::Invalid),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Key::Id.to_string(), "ID");
        assert_eq!(Key::Number.to_string(), "Number");
        assert_eq!(Key::Type.to_string(), "Type");
        assert_eq!(Key::Description.to_string(), "Description");
        assert_eq!(Key::Idx.to_string(), "IDX");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("ID".parse(), Ok(Key::Id));
        assert_eq!("Number".parse(), Ok(Key::Number));
        assert_eq!("Type".parse(), Ok(Key::Type));
        assert_eq!("Description".parse(), Ok(Key::Description));
        assert_eq!("IDX".parse(), Ok(Key::Idx));

        assert_eq!("".parse::<Key>(), Err(ParseError::Empty));
        assert_eq!("Noodles".parse::<Key>(), Err(ParseError::Invalid));
    }
}
