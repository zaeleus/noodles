use std::{error, fmt, mem};

use indexmap::IndexSet;

use crate::variant::record_buf::Filters;

/// An error when raw VCF record filters fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A filter is invalid.
    InvalidFilter,
    /// A filter is duplicated.
    DuplicateFilter,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::InvalidFilter => write!(f, "invalid filter"),
            Self::DuplicateFilter => write!(f, "duplicate filter"),
        }
    }
}

pub(super) fn parse_filters(s: &str, filters: &mut Option<Filters>) -> Result<(), ParseError> {
    const DELIMITER: char = ';';
    const PASS: &str = "PASS";

    if s.is_empty() {
        return Err(ParseError::Empty);
    } else if s == PASS {
        *filters = Some(Filters::pass());
        return Ok(());
    }

    let mut set = match mem::take(filters) {
        None => IndexSet::new(),
        Some(filters) => {
            let mut set: IndexSet<String> = filters.into();
            set.clear();
            set
        }
    };

    for raw_filter in s.split(DELIMITER) {
        if !set.insert(raw_filter.into()) {
            return Err(ParseError::DuplicateFilter);
        } else if !is_valid_filter(raw_filter) {
            return Err(ParseError::InvalidFilter);
        }
    }

    *filters = Some(Filters(set));

    Ok(())
}

fn is_valid_filter(s: &str) -> bool {
    match s {
        "" | "0" => false,
        _ => s.chars().all(|c| !c.is_whitespace()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_filters() -> Result<(), ParseError> {
        let mut filters = None;

        parse_filters("PASS", &mut filters)?;
        assert_eq!(filters, Some(Filters::pass()));

        parse_filters("q10", &mut filters)?;
        assert_eq!(filters, Some([String::from("q10")].into_iter().collect()));

        parse_filters("q10;s50", &mut filters)?;
        assert_eq!(
            filters,
            Some(
                [String::from("q10"), String::from("s50")]
                    .into_iter()
                    .collect()
            )
        );

        assert_eq!(parse_filters("", &mut filters), Err(ParseError::Empty));
        assert_eq!(
            parse_filters("0", &mut filters),
            Err(ParseError::InvalidFilter)
        );
        assert_eq!(
            parse_filters("q 10", &mut filters),
            Err(ParseError::InvalidFilter)
        );
        assert_eq!(
            parse_filters(";q10", &mut filters),
            Err(ParseError::InvalidFilter)
        );
        assert_eq!(
            parse_filters("q10;;s50", &mut filters),
            Err(ParseError::InvalidFilter)
        );
        assert_eq!(
            parse_filters("q10;", &mut filters),
            Err(ParseError::InvalidFilter)
        );
        assert_eq!(
            parse_filters("q10;q10", &mut filters),
            Err(ParseError::DuplicateFilter)
        );

        Ok(())
    }
}
