//! GTF lines.

use std::{error, fmt, io, str::FromStr};

use super::{record_buf, Line, RecordBuf};

const COMMENT_PREFIX: char = '#';

/// A GTF line.
#[derive(Clone, Debug, PartialEq)]
pub enum LineBuf {
    /// A comment (`#`).
    Comment(String),
    /// A record.
    Record(RecordBuf),
}

impl TryFrom<Line> for LineBuf {
    type Error = io::Error;

    fn try_from(line: Line) -> Result<Self, Self::Error> {
        if let Some(s) = line.as_comment() {
            Ok(Self::Comment(s.into()))
        } else if let Some(result) = line.as_record() {
            result.and_then(RecordBuf::try_from).map(Self::Record)
        } else {
            unreachable!()
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
    fn test_from_str() {
        assert_eq!(
            "##format: gtf".parse(),
            Ok(LineBuf::Comment(String::from("#format: gtf")))
        );

        let s = "sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id \"ndls0\"; transcript_id \"ndls0\";";
        assert!(matches!(s.parse(), Ok(LineBuf::Record(_))));
    }
}
