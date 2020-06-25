//! GFF lines.

use std::{error, fmt, str::FromStr};

use super::{directive, record, Directive, Record};

/// A GFF line.
#[derive(Clone, Debug, PartialEq)]
pub enum Line {
    /// A directive (`##`).
    Directive(Directive),
    /// A comment (`#`),
    Comment(String),
    /// A record.
    Record(Record),
}

/// An error returns when a raw GFF line fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The directive is invalid.
    InvalidDirective(directive::ParseError),
    /// The record is invalid.
    InvalidRecord(record::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDirective(e) => write!(f, "{}", e),
            Self::InvalidRecord(e) => write!(f, "{}", e),
        }
    }
}

impl FromStr for Line {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.starts_with(directive::PREFIX) {
            s.parse()
                .map(Self::Directive)
                .map_err(ParseError::InvalidDirective)
        } else if s.starts_with('#') {
            Ok(Self::Comment(s[1..].into()))
        } else {
            s.parse()
                .map(Self::Record)
                .map_err(ParseError::InvalidRecord)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let s = "##gff-version 3";
        let actual = s.parse::<Line>()?;
        let expected = Line::Directive(Directive::GffVersion(String::from("3")));
        assert_eq!(actual, expected);

        let s = "#format: gff3";
        let actual = s.parse::<Line>()?;
        let expected = Line::Comment(String::from("format: gff3"));
        assert_eq!(actual, expected);

        let s = "sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id=ndls0;gene_name=gene0";
        let actual = s.parse::<Line>()?;
        assert!(matches!(actual, Line::Record(_)));

        Ok(())
    }
}
