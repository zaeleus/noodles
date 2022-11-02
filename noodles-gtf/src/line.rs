//! GTF lines.

use std::{error, fmt, str::FromStr};

use super::{record, Record};

/// A GTF line.
#[derive(Clone, Debug, PartialEq)]
pub enum Line {
    /// A comment (`#`).
    Comment(String),
    /// A record.
    Record(Record),
}

impl fmt::Display for Line {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Line::Comment(comment) => write!(f, "#{}", comment),
            Line::Record(record) => write!(f, "{}", record),
        }
    }
}

/// An error returns when a raw GFF line fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The record is invalid.
    InvalidRecord(record::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord(e) => write!(f, "invalid record: {}", e),
        }
    }
}

impl FromStr for Line {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some(t) = s.strip_prefix('#') {
            Ok(Self::Comment(t.into()))
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
    fn test_fmt() {
        let line = Line::Comment(String::from("format: gtf"));
        assert_eq!(line.to_string(), "#format: gtf");

        let line = Line::Record(Record::default());
        assert_eq!(line.to_string(), ".\t.\t.\t1\t1\t.\t.\t.\t");
    }

    #[test]
    fn test_from_str() {
        assert_eq!(
            "##format: gtf".parse(),
            Ok(Line::Comment(String::from("#format: gtf")))
        );

        let s = "sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id \"ndls0\"; transcript_id \"ndls0\";";
        assert!(matches!(s.parse(), Ok(Line::Record(_))));
    }
}
