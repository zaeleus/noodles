mod kind;

use std::{error, fmt, str::FromStr};

pub use self::kind::Kind;

type Field = (String, String);

const DELIMITER: char = '\t';
const DATA_FIELD_DELIMITER: char = ':';

/// A SAM header record.
#[derive(Debug, Eq, PartialEq)]
pub enum Record {
    Header(Vec<Field>),
    ReferenceSequence(Vec<Field>),
    ReadGroup(Vec<Field>),
    Program(Vec<Field>),
    Comment(String),
}

/// An error returned when a raw SAM header record fails to parse.
#[derive(Debug)]
pub enum ParseError {
    /// The kind is missing.
    MissingKind,
    /// The kind is invalid.
    InvalidKind(kind::ParseError),
    /// A tag is missing.
    MissingTag,
    /// A tag value is missing.
    MissingValue(String),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid record: ")?;

        match self {
            Self::MissingKind => write!(f, "missing kind"),
            Self::InvalidKind(e) => write!(f, "{}", e),
            Self::MissingTag => write!(f, "missing field tag"),
            Self::MissingValue(tag) => write!(f, "missing value for {}", tag),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut pieces = s.split(DELIMITER);

        let kind = pieces
            .next()
            .ok_or_else(|| ParseError::MissingKind)
            .and_then(|s| s.parse().map_err(ParseError::InvalidKind))?;

        let record = if let Kind::Comment = kind {
            pieces
                .next()
                .map(|s| Self::Comment(s.into()))
                .ok_or_else(|| ParseError::MissingValue(Kind::Comment.to_string()))?
        } else {
            let fields = pieces
                .map(|field| {
                    let mut field_pieces = field.splitn(2, DATA_FIELD_DELIMITER);

                    let tag = field_pieces.next().ok_or_else(|| ParseError::MissingTag)?;
                    let value = field_pieces
                        .next()
                        .ok_or_else(|| ParseError::MissingValue(tag.into()))?;

                    Ok((tag.into(), value.into()))
                })
                .collect::<Result<_, _>>()?;

            match kind {
                Kind::Header => Self::Header(fields),
                Kind::ReferenceSequence => Self::ReferenceSequence(fields),
                Kind::ReadGroup => Self::ReadGroup(fields),
                Kind::Program => Self::Program(fields),
                Kind::Comment => unreachable!(),
            }
        };

        Ok(record)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let record = "@HD\tVN:1.6".parse()?;
        let fields = vec![(String::from("VN"), String::from("1.6"))];
        assert!(matches!(record, Record::Header(f) if f == fields));

        let record = "@SQ\tSN:sq0\tLN:8".parse()?;
        let fields = vec![
            (String::from("SN"), String::from("sq0")),
            (String::from("LN"), String::from("8")),
        ];
        assert!(matches!(record, Record::ReferenceSequence(f) if f == fields));

        let record = "@RG\tID:rg0".parse()?;
        let fields = vec![(String::from("ID"), String::from("rg0"))];
        assert!(matches!(record, Record::ReadGroup(f) if f == fields));

        let record = "@PG\tID:pg0".parse()?;
        let fields = vec![(String::from("ID"), String::from("pg0"))];
        assert!(matches!(record, Record::Program(f) if f == fields));

        let record = "@CO\tnoodles".parse()?;
        assert!(matches!(record, Record::Comment(c) if c == "noodles"));

        Ok(())
    }

    #[test]
    fn test_from_str_when_comment_has_no_message() {
        assert!("@CO".parse::<Record>().is_err());
    }
}
