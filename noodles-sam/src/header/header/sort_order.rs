use std::{error, fmt, str::FromStr};

/// A SAM header header sort order (`SO`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SortOrder {
    /// Alignment order is unknown (`unknown`).
    Unknown,
    /// Alignments are not sorted (`unsorted`).
    Unsorted,
    /// Alignments are sorted by read name (`queryname`).
    QueryName,
    /// Alignments are sorted by reference sequence and position (`coordinate`).
    Coordinate,
}

impl AsRef<str> for SortOrder {
    fn as_ref(&self) -> &str {
        match self {
            Self::Unknown => "unknown",
            Self::Unsorted => "unsorted",
            Self::QueryName => "queryname",
            Self::Coordinate => "coordinate",
        }
    }
}

impl Default for SortOrder {
    fn default() -> Self {
        Self::Unknown
    }
}

impl fmt::Display for SortOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

/// An error returned when a raw SAM header header sort order fails to parse.
#[derive(Debug)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "invalid sort order: expected {{unknown, unsorted, queryname, coordinate}}, got {}",
            self.0
        )
    }
}

impl FromStr for SortOrder {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "unknown" => Ok(Self::Unknown),
            "unsorted" => Ok(Self::Unsorted),
            "queryname" => Ok(Self::QueryName),
            "coordinate" => Ok(Self::Coordinate),
            _ => Err(ParseError(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(SortOrder::default(), SortOrder::Unknown);
    }

    #[test]
    fn test_fmt() {
        assert_eq!(SortOrder::Unknown.to_string(), "unknown");
        assert_eq!(SortOrder::Unsorted.to_string(), "unsorted");
        assert_eq!(SortOrder::QueryName.to_string(), "queryname");
        assert_eq!(SortOrder::Coordinate.to_string(), "coordinate");
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("unknown".parse::<SortOrder>()?, SortOrder::Unknown);
        assert_eq!("unsorted".parse::<SortOrder>()?, SortOrder::Unsorted);
        assert_eq!("queryname".parse::<SortOrder>()?, SortOrder::QueryName);
        assert_eq!("coordinate".parse::<SortOrder>()?, SortOrder::Coordinate);

        assert!("".parse::<SortOrder>().is_err());
        assert!("noodles".parse::<SortOrder>().is_err());
        assert!("QueryName".parse::<SortOrder>().is_err());

        Ok(())
    }
}
