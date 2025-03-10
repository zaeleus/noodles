//! GTF lines.

use std::{error, fmt, str::FromStr};

use super::{record_buf, RecordBuf};

const COMMENT_PREFIX: char = '#';

/// A GTF line.
#[derive(Clone, Debug, PartialEq)]
pub enum LineBuf {
    /// A comment (`#`).
    Comment(String),
    /// A record.
    Record(RecordBuf),
}

impl fmt::Display for LineBuf {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Comment(comment) => write!(f, "{COMMENT_PREFIX}{comment}"),
            Self::Record(record) => write!(f, "{record}"),
        }
    }
}

/// An error returns when a raw GFF line fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The record is invalid.
    InvalidRecord(record_buf::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidRecord(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord(_) => f.write_str("invalid record"),
        }
    }
}

impl FromStr for LineBuf {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some(t) = s.strip_prefix(COMMENT_PREFIX) {
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
        let line = LineBuf::Comment(String::from("format: gtf"));
        assert_eq!(line.to_string(), "#format: gtf");

        let line = LineBuf::Record(RecordBuf::default());
        assert_eq!(line.to_string(), ".\t.\t.\t1\t1\t.\t.\t.\t");
    }

    #[test]
    fn test_from_str() {
        assert_eq!(
            "##format: gtf".parse(),
            Ok(LineBuf::Comment(String::from("#format: gtf")))
        );

        let s = "sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id \"ndls0\"; transcript_id \"ndls0\";";
        assert!(matches!(s.parse(), Ok(LineBuf::Record(_))));
    }
}
