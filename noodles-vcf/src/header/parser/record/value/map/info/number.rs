use std::{error, fmt, num};

use crate::header::Number;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid(num::ParseIntError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Empty => None,
            Self::Invalid(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid(_) => f.write_str("invalid input"),
        }
    }
}

pub(super) fn parse_number(s: &str) -> Result<Number, ParseError> {
    match s {
        "" => Err(ParseError::Empty),
        "A" => Ok(Number::AlternateBases),
        "R" => Ok(Number::ReferenceAlternateBases),
        "G" => Ok(Number::Samples),
        "." => Ok(Number::Unknown),
        _ => s.parse().map(Number::Count).map_err(ParseError::Invalid),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_number() {
        assert_eq!(parse_number("1"), Ok(Number::Count(1)));
        assert_eq!(parse_number("A"), Ok(Number::AlternateBases));
        assert_eq!(parse_number("R"), Ok(Number::ReferenceAlternateBases));
        assert_eq!(parse_number("G"), Ok(Number::Samples));
        assert_eq!(parse_number("."), Ok(Number::Unknown));

        assert_eq!(parse_number(""), Err(ParseError::Empty));
        assert!(matches!(parse_number("ndls"), Err(ParseError::Invalid(_))));
    }
}
