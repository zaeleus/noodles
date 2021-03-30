//! SAM header header subsort order.

use std::{error, fmt, str::FromStr};

use super::SortOrder;

const DELIMITER: char = ':';

/// A SAM header header subsort order (`SS`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SubsortOrder {
    /// Alignments are primarily unsorted (`unsorted`).
    Unsorted(Vec<String>),
    /// Alignments are primarily sorted by read name (`queryname`).
    QueryName(Vec<String>),
    /// Alignments are primarily sorted by reference sequence and position (`coordinate`).
    Coordinate(Vec<String>),
}

impl fmt::Display for SubsortOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unsorted(subsorts) => write_orders(f, SortOrder::Unsorted, subsorts),
            Self::QueryName(subsorts) => write_orders(f, SortOrder::QueryName, subsorts),
            Self::Coordinate(subsorts) => write_orders(f, SortOrder::Coordinate, subsorts),
        }
    }
}

fn write_orders(f: &mut fmt::Formatter<'_>, sort: SortOrder, subsorts: &[String]) -> fmt::Result {
    write!(f, "{}", sort)?;

    for subsort in subsorts {
        write!(f, "{}{}", DELIMITER, subsort)?;
    }

    Ok(())
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

        let mut pieces = s.splitn(2, DELIMITER);

        let order = pieces.next().ok_or(ParseError::MissingOrder)?;

        let subsort = pieces
            .next()
            .map(|s| s.split(DELIMITER).map(|t| t.into()).collect())
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
            SubsortOrder::Unsorted(vec![String::from("MI")]).to_string(),
            "unsorted:MI"
        );

        assert_eq!(
            SubsortOrder::QueryName(vec![String::from("MI")]).to_string(),
            "queryname:MI"
        );

        assert_eq!(
            SubsortOrder::Coordinate(vec![String::from("MI")]).to_string(),
            "coordinate:MI"
        );

        assert_eq!(
            SubsortOrder::Unsorted(vec![String::from("MI"), String::from("coordinate")])
                .to_string(),
            "unsorted:MI:coordinate"
        );
    }

    #[test]
    fn test_from_str() {
        assert_eq!(
            "unsorted:MI".parse(),
            Ok(SubsortOrder::Unsorted(vec![String::from("MI")]))
        );
        assert_eq!(
            "queryname:MI".parse(),
            Ok(SubsortOrder::QueryName(vec![String::from("MI")]))
        );
        assert_eq!(
            "coordinate:MI".parse(),
            Ok(SubsortOrder::Coordinate(vec![String::from("MI")]))
        );
        assert_eq!(
            "unsorted:MI:coordinate".parse(),
            Ok(SubsortOrder::Unsorted(vec![
                String::from("MI"),
                String::from("coordinate")
            ]))
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
