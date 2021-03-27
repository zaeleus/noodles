//! SAM header header group order.

use std::{error, fmt, str::FromStr};

/// A SAM header header group order (`GO`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GroupOrder {
    /// Alignments are not grouped (`none`).
    None,
    /// Alignments are grouped by read name (`query`).
    Query,
    /// Alignments are grouped by reference sequence and position (`reference`).
    Reference,
}

impl AsRef<str> for GroupOrder {
    fn as_ref(&self) -> &str {
        match self {
            Self::None => "none",
            Self::Query => "query",
            Self::Reference => "reference",
        }
    }
}

impl Default for GroupOrder {
    fn default() -> Self {
        Self::None
    }
}

impl fmt::Display for GroupOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

/// An error returned when a raw SAM header header group order fails to parse.
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

impl FromStr for GroupOrder {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            "none" => Ok(Self::None),
            "query" => Ok(Self::Query),
            "reference" => Ok(Self::Reference),
            _ => Err(ParseError::Invalid),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(GroupOrder::default(), GroupOrder::None);
    }

    #[test]
    fn test_fmt() {
        assert_eq!(GroupOrder::None.to_string(), "none");
        assert_eq!(GroupOrder::Query.to_string(), "query");
        assert_eq!(GroupOrder::Reference.to_string(), "reference");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("none".parse(), Ok(GroupOrder::None));
        assert_eq!("query".parse(), Ok(GroupOrder::Query));
        assert_eq!("reference".parse(), Ok(GroupOrder::Reference));

        assert_eq!("".parse::<GroupOrder>(), Err(ParseError::Empty));
        assert_eq!("noodles".parse::<GroupOrder>(), Err(ParseError::Invalid));
        assert_eq!("Query".parse::<GroupOrder>(), Err(ParseError::Invalid));
    }
}
