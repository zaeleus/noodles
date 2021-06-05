//! SAM header record and components.

pub mod kind;
pub mod value;

use std::{error, fmt, str::FromStr};

use indexmap::IndexMap;

pub use self::{kind::Kind, value::Value};

const DELIMITER: char = '\t';
const DATA_FIELD_DELIMITER: char = ':';

/// A generic SAM header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Record {
    kind: Kind,
    value: Value,
}

impl Record {
    /// Creates a generic SAM header record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{record::{Kind, Value}, Record};
    /// let record = Record::new(Kind::Comment, Value::String(String::from("noodles-sam")));
    /// ```
    pub fn new(kind: Kind, value: Value) -> Self {
        Self { kind, value }
    }

    /// Returns the kind of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{record::{Kind, Value}, Record};
    /// let record = Record::new(Kind::Comment, Value::String(String::from("noodles-sam")));
    /// assert_eq!(record.kind(), Kind::Comment);
    /// ```
    pub fn kind(&self) -> Kind {
        self.kind
    }

    /// Returns the value of the record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{record::{Kind, Value}, Record};
    /// let record = Record::new(Kind::Comment, Value::String(String::from("noodles-sam")));
    /// assert_eq!(record.value(), &Value::String(String::from("noodles-sam")));
    /// ```
    pub fn value(&self) -> &Value {
        &self.value
    }
}

impl From<Record> for (Kind, Value) {
    fn from(record: Record) -> Self {
        (record.kind, record.value)
    }
}

/// An error returned when a raw SAM header record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The kind is missing.
    MissingKind,
    /// The kind is invalid.
    InvalidKind(kind::ParseError),
    /// A tag is missing.
    MissingTag,
    /// A tag is duplicated.
    DuplicateTag(String),
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
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {}", tag),
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
            .ok_or(ParseError::MissingKind)
            .and_then(|s| s.parse().map_err(ParseError::InvalidKind))?;

        let value = match kind {
            Kind::Comment => parse_comment(&mut pieces)?,
            _ => parse_map(&mut pieces)?,
        };

        Ok(Self::new(kind, value))
    }
}

fn parse_comment<'a, I>(iter: &mut I) -> Result<Value, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    iter.next()
        .map(|s| Value::String(s.into()))
        .ok_or_else(|| ParseError::MissingValue(Kind::Comment.to_string()))
}

fn parse_map<'a, I>(iter: &mut I) -> Result<Value, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    let mut map = IndexMap::new();

    for s in iter {
        let mut components = s.splitn(2, DATA_FIELD_DELIMITER);

        let tag = components.next().ok_or(ParseError::MissingTag)?;
        let value = components
            .next()
            .ok_or_else(|| ParseError::MissingValue(tag.into()))?;

        if map.insert(tag.into(), value.into()).is_some() {
            return Err(ParseError::DuplicateTag(tag.into()));
        }
    }

    Ok(Value::Map(map))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), value::TryFromIteratorError> {
        assert_eq!(
            "@HD\tVN:1.6".parse(),
            Ok(Record::new(
                Kind::Header,
                Value::try_from_iter(vec![("VN", "1.6")])?,
            ))
        );

        assert_eq!(
            "@SQ\tSN:sq0\tLN:8".parse(),
            Ok(Record::new(
                Kind::ReferenceSequence,
                Value::try_from_iter(vec![("SN", "sq0"), ("LN", "8"),])?
            ))
        );

        assert_eq!(
            "@RG\tID:rg0".parse(),
            Ok(Record::new(
                Kind::ReadGroup,
                Value::try_from_iter(vec![("ID", "rg0")])?
            ))
        );

        assert_eq!(
            "@PG\tID:pg0".parse(),
            Ok(Record::new(
                Kind::Program,
                Value::try_from_iter(vec![("ID", "pg0")])?
            ))
        );

        assert_eq!(
            "@CO\tnoodles-sam".parse(),
            Ok(Record::new(
                Kind::Comment,
                Value::String(String::from("noodles-sam"))
            ))
        );

        assert_eq!(
            "@CO\t".parse(),
            Ok(Record::new(Kind::Comment, Value::String(String::from(""))))
        );

        assert_eq!(
            "@CO".parse::<Record>(),
            Err(ParseError::MissingValue(String::from("@CO")))
        );

        Ok(())
    }
}
