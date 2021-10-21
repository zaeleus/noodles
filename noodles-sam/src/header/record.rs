//! SAM header record and components.

pub mod kind;
pub mod value;

use std::{error, fmt, str::FromStr};

pub use self::{kind::Kind, value::Value};

use self::value::Fields;

const DELIMITER: char = '\t';
const DATA_FIELD_DELIMITER: char = ':';
const TAG_LENGTH: usize = 2;

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
    /// The input is invalid.
    Invalid,
    /// The kind is missing.
    MissingKind,
    /// The kind is invalid.
    InvalidKind(kind::ParseError),
    /// A field is invalid.
    InvalidField,
    /// A tag is invalid.
    InvalidTag,
    /// A value is invalid.
    InvalidValue,
    /// A tag is duplicated.
    DuplicateTag(String),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => write!(f, "invalid input"),
            Self::MissingKind => write!(f, "missing kind"),
            Self::InvalidKind(e) => write!(f, "invalid kind: {}", e),
            Self::InvalidField => write!(f, "invalid field"),
            Self::InvalidTag => write!(f, "invalid tag"),
            Self::InvalidValue => write!(f, "invalid value"),
            Self::DuplicateTag(tag) => write!(f, "duplicate tag: {}", tag),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.split_once(DELIMITER) {
            Some((k, v)) => {
                let kind = k.parse().map_err(ParseError::InvalidKind)?;

                let value = match kind {
                    Kind::Comment => parse_comment(v)?,
                    _ => parse_map(v)?,
                };

                Ok(Self::new(kind, value))
            }
            None => Err(ParseError::Invalid),
        }
    }
}

fn parse_comment(s: &str) -> Result<Value, ParseError> {
    Ok(Value::String(s.into()))
}

fn parse_map(s: &str) -> Result<Value, ParseError> {
    let raw_fields = s.split(DELIMITER);
    let mut fields = Fields::new();

    for raw_field in raw_fields {
        let (tag, value) = parse_field(raw_field)?;

        if fields.insert(tag.clone(), value).is_some() {
            return Err(ParseError::DuplicateTag(tag));
        }
    }

    Ok(Value::Map(fields))
}

fn parse_field(s: &str) -> Result<(String, String), ParseError> {
    match s.split_once(DATA_FIELD_DELIMITER) {
        Some((tag, value)) => {
            if !is_valid_tag(tag) {
                return Err(ParseError::InvalidTag);
            }

            if !is_valid_value(value) {
                return Err(ParseError::InvalidValue);
            }

            Ok((tag.into(), value.into()))
        }
        None => Err(ParseError::InvalidField),
    }
}

fn is_valid_tag(s: &str) -> bool {
    if s.len() != TAG_LENGTH {
        return false;
    }

    let mut chars = s.chars();

    if let Some(c) = chars.next() {
        if !c.is_ascii_alphabetic() {
            return false;
        }
    }

    if let Some(c) = chars.next() {
        if !c.is_ascii_alphanumeric() {
            return false;
        }
    }

    true
}

fn is_valid_value(s: &str) -> bool {
    !s.is_empty() && s.chars().all(|c| matches!(c, ' '..='~'))
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
                Value::try_from_iter([("VN", "1.6")])?,
            ))
        );

        assert_eq!(
            "@SQ\tSN:sq0\tLN:8".parse(),
            Ok(Record::new(
                Kind::ReferenceSequence,
                Value::try_from_iter([("SN", "sq0"), ("LN", "8"),])?
            ))
        );

        assert_eq!(
            "@RG\tID:rg0".parse(),
            Ok(Record::new(
                Kind::ReadGroup,
                Value::try_from_iter([("ID", "rg0")])?
            ))
        );

        assert_eq!(
            "@PG\tID:pg0".parse(),
            Ok(Record::new(
                Kind::Program,
                Value::try_from_iter([("ID", "pg0")])?
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

        assert!(matches!(
            "@ND\t".parse::<Record>(),
            Err(ParseError::InvalidKind(_))
        ));

        assert_eq!("@HD\t".parse::<Record>(), Err(ParseError::InvalidField));
        assert_eq!("@HD\tVN".parse::<Record>(), Err(ParseError::InvalidField));

        assert_eq!("@HD\tV:1.6".parse::<Record>(), Err(ParseError::InvalidTag));
        assert_eq!("@HD\t0V:1.6".parse::<Record>(), Err(ParseError::InvalidTag));
        assert_eq!(
            "@HD\tVER:1.6".parse::<Record>(),
            Err(ParseError::InvalidTag)
        );

        assert_eq!("@PG\tID:".parse::<Record>(), Err(ParseError::InvalidValue));
        assert_eq!(
            "@PG\tID:🍜".parse::<Record>(),
            Err(ParseError::InvalidValue)
        );

        assert_eq!(
            "@HD\tVN:1.6\tVN:1.6".parse::<Record>(),
            Err(ParseError::DuplicateTag(String::from("VN")))
        );

        assert_eq!("@CO".parse::<Record>(), Err(ParseError::Invalid));

        Ok(())
    }
}
