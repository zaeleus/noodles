//! SAM header header subsort order.

use std::{error, fmt, str::FromStr};

/// A SAM header header subsort order (`SS`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SubsortOrder {
    /// Alignments are primarily unsorted (`unsorted`).
    Unsorted(String),
    /// Alignments are primarily sorted by read name (`queryname`).
    QueryName(String),
    /// Alignments are primarily sorted by reference sequence and position (`coordinate`).
    Coordinate(String),
}

impl fmt::Display for SubsortOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unsorted(subsort) => write!(f, "unsorted:{}", subsort),
            Self::QueryName(subsort) => write!(f, "queryname:{}", subsort),
            Self::Coordinate(subsort) => write!(f, "coordinate:{}", subsort),
        }
    }
}

/// An error returned when a raw SAM header header subsort order fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The primary sort order is missing.
    MissingOrder,
    /// The primary sort order is invalid.
    InvalidOrder,
    /// The subsort order is missing.
    MissingSubsort,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::MissingOrder => f.write_str("missing primary sort order"),
            Self::InvalidOrder => f.write_str("invalid primary sort order"),
            Self::MissingSubsort => f.write_str("missing subsort order"),
        }
    }
}

impl FromStr for SubsortOrder {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let mut pieces = s.splitn(2, ':');

        let order = pieces.next().ok_or(ParseError::MissingOrder)?;

        let subsort = pieces
            .next()
            .map(|s| s.into())
            .ok_or(ParseError::MissingSubsort)?;

        match order {
            "unsorted" => Ok(Self::Unsorted(subsort)),
            "queryname" => Ok(Self::QueryName(subsort)),
            "coordinate" => Ok(Self::Coordinate(subsort)),
            _ => Err(ParseError::InvalidOrder),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(
            SubsortOrder::Unsorted(String::from("MI")).to_string(),
            "unsorted:MI"
        );

        assert_eq!(
            SubsortOrder::QueryName(String::from("MI")).to_string(),
            "queryname:MI"
        );

        assert_eq!(
            SubsortOrder::Coordinate(String::from("MI")).to_string(),
            "coordinate:MI"
        );
    }

    #[test]
    fn test_from_str() {
        assert_eq!(
            "unsorted:MI".parse(),
            Ok(SubsortOrder::Unsorted(String::from("MI")))
        );
        assert_eq!(
            "queryname:MI".parse(),
            Ok(SubsortOrder::QueryName(String::from("MI")))
        );
        assert_eq!(
            "coordinate:MI".parse(),
            Ok(SubsortOrder::Coordinate(String::from("MI")))
        );
        assert_eq!(
            "unsorted:MI:coordinate".parse(),
            Ok(SubsortOrder::Unsorted(String::from("MI:coordinate")))
        );

        assert_eq!("".parse::<SubsortOrder>(), Err(ParseError::Empty));
        assert_eq!(
            "noodles".parse::<SubsortOrder>(),
            Err(ParseError::MissingSubsort)
        );
        assert_eq!(
            "QueryName".parse::<SubsortOrder>(),
            Err(ParseError::MissingSubsort)
        );
        assert_eq!(
            "queryname".parse::<SubsortOrder>(),
            Err(ParseError::MissingSubsort)
        );
    }
}
