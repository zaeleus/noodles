//! VCF record filters.

use std::{error, fmt, str::FromStr};

use indexmap::IndexSet;

use super::MISSING_FIELD;

const PASS_STATUS: &str = "PASS";
const DELIMITER: char = ';';

/// VCF record filters (`FILTER`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Filters {
    /// Missing (`*`).
    Missing,
    /// Pass (`PASS`).
    Pass,
    /// A list of filters that caused the record to fail.
    Fail(IndexSet<String>),
}

impl Default for Filters {
    fn default() -> Self {
        Self::Missing
    }
}

impl fmt::Display for Filters {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Missing => f.write_str(MISSING_FIELD),
            Self::Pass => f.write_str(PASS_STATUS),
            Self::Fail(ids) => {
                for (i, id) in ids.iter().enumerate() {
                    if i > 0 {
                        write!(f, "{}", DELIMITER)?;
                    }

                    f.write_str(id)?;
                }

                Ok(())
            }
        }
    }
}

/// An error returned when a raw VCF filter fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A filter is duplicated.
    DuplicateFilter(String),
    /// A filter is invalid.
    InvalidFilter(String),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::DuplicateFilter(filter) => write!(f, "duplicate filter: {}", filter),
            Self::InvalidFilter(s) => write!(f, "invalid filter: {}", s),
        }
    }
}

impl FromStr for Filters {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "" => Err(ParseError::Empty),
            MISSING_FIELD => Ok(Self::Missing),
            PASS_STATUS => Ok(Self::Pass),
            _ => {
                let mut filters = IndexSet::new();

                for filter in s.split(DELIMITER) {
                    if !filters.insert(filter.into()) {
                        return Err(ParseError::DuplicateFilter(filter.into()));
                    } else if !is_valid_filter(filter) {
                        return Err(ParseError::InvalidFilter(filter.into()));
                    }
                }

                Ok(Self::Fail(filters))
            }
        }
    }
}

fn is_valid_filter(s: &str) -> bool {
    match s {
        "" | "0" => false,
        _ => s.chars().all(|c| !c.is_ascii_whitespace()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(Filters::default(), Filters::Missing);
    }

    #[test]
    fn test_fmt() {
        assert_eq!(Filters::Missing.to_string(), ".");
        assert_eq!(Filters::Pass.to_string(), "PASS");

        let filters = Filters::Fail(vec![String::from("q10")].into_iter().collect());
        assert_eq!(filters.to_string(), "q10");

        let filters = Filters::Fail(
            vec![String::from("q10"), String::from("s50")]
                .into_iter()
                .collect(),
        );
        assert_eq!(filters.to_string(), "q10;s50");
    }

    #[test]
    fn test_from_str() {
        assert_eq!(".".parse(), Ok(Filters::Missing));
        assert_eq!("PASS".parse(), Ok(Filters::Pass));

        assert_eq!(
            "q10".parse(),
            Ok(Filters::Fail(
                vec![String::from("q10")].into_iter().collect()
            ))
        );

        assert_eq!(
            "q10;s50".parse(),
            Ok(Filters::Fail(
                vec![String::from("q10"), String::from("s50")]
                    .into_iter()
                    .collect()
            ))
        );

        assert_eq!("".parse::<Filters>(), Err(ParseError::Empty));
        assert_eq!(
            "q10;q10".parse::<Filters>(),
            Err(ParseError::DuplicateFilter(String::from("q10")))
        );
        assert_eq!(
            "0".parse::<Filters>(),
            Err(ParseError::InvalidFilter(String::from("0")))
        );
        assert_eq!(
            "q 10".parse::<Filters>(),
            Err(ParseError::InvalidFilter(String::from("q 10")))
        );
        assert_eq!(
            ";q10".parse::<Filters>(),
            Err(ParseError::InvalidFilter(String::from("")))
        );
        assert_eq!(
            "q10;;s50".parse::<Filters>(),
            Err(ParseError::InvalidFilter(String::from("")))
        );
        assert_eq!(
            "q10;".parse::<Filters>(),
            Err(ParseError::InvalidFilter(String::from("")))
        );
    }
}
