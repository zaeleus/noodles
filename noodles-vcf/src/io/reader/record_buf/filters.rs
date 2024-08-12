use std::{error, fmt};

use crate::variant::record_buf::Filters;

/// An error when raw VCF record filters fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A filter is duplicated.
    DuplicateFilter,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::DuplicateFilter => write!(f, "duplicate filter"),
        }
    }
}

pub(super) fn parse_filters(s: &str, filters: &mut Filters) -> Result<(), ParseError> {
    const DELIMITER: char = ';';
    const PASS: &str = "PASS";

    if s.is_empty() {
        return Err(ParseError::Empty);
    } else if s == PASS {
        *filters = Filters::pass();
        return Ok(());
    }

    let filters = filters.as_mut();
    filters.clear();

    for raw_filter in s.split(DELIMITER) {
        if !filters.insert(raw_filter.into()) {
            return Err(ParseError::DuplicateFilter);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_filters() -> Result<(), ParseError> {
        let mut filters = Filters::default();

        parse_filters("PASS", &mut filters)?;
        assert_eq!(filters, Filters::pass());

        parse_filters("q10", &mut filters)?;
        assert_eq!(filters, [String::from("q10")].into_iter().collect());

        parse_filters("q10;s50", &mut filters)?;
        assert_eq!(
            filters,
            [String::from("q10"), String::from("s50")]
                .into_iter()
                .collect()
        );

        assert_eq!(parse_filters("", &mut filters), Err(ParseError::Empty));

        assert_eq!(
            parse_filters("q10;q10", &mut filters),
            Err(ParseError::DuplicateFilter)
        );

        Ok(())
    }
}
