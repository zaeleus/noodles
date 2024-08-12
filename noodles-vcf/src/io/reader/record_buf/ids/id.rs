use std::{error, fmt};

/// An error when a raw VCF record ID fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
        }
    }
}

pub(super) fn parse_id(s: &str) -> Result<&str, ParseError> {
    if s.is_empty() {
        Err(ParseError::Empty)
    } else {
        Ok(s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_id() {
        assert_eq!(parse_id("nd0"), Ok("nd0"));
        assert_eq!(parse_id(""), Err(ParseError::Empty));
    }
}
