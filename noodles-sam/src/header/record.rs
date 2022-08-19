//! SAM header record.

pub mod kind;
pub mod value;

use std::{error, fmt, str::FromStr};

pub use self::kind::Kind;

use self::value::{
    map::{self, Program, ReadGroup, ReferenceSequence},
    Map,
};

const DELIMITER: char = '\t';
const TAG_LENGTH: usize = 2;

/// A SAM header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Record {
    /// A header (`HD`) record.
    Header(Map<map::Header>),
    /// A reference sequence (`SQ`) record.
    ReferenceSequence(Map<ReferenceSequence>),
    /// A read group (`RG`) record.
    ReadGroup(Map<ReadGroup>),
    /// A program (`PG`) record.
    Program(Map<Program>),
    /// A comment (`CO`) record.
    Comment(String),
}

/// An error returned when a raw SAM header record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
    /// The kind is invalid.
    InvalidKind(kind::ParseError),
    /// A field is invalid.
    InvalidField,
    /// A tag is invalid.
    InvalidTag,
    /// A value is invalid.
    InvalidValue,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => write!(f, "invalid input"),
            Self::InvalidKind(e) => write!(f, "invalid kind: {}", e),
            Self::InvalidField => write!(f, "invalid field"),
            Self::InvalidTag => write!(f, "invalid tag"),
            Self::InvalidValue => write!(f, "invalid value"),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (k, v) = s.split_once(DELIMITER).ok_or(ParseError::Invalid)?;
        let kind = k.parse().map_err(ParseError::InvalidKind)?;

        match kind {
            Kind::Header => {
                let fields = split_fields(v)?;
                let header =
                    Map::<map::Header>::try_from(fields).map_err(|_| ParseError::Invalid)?;
                Ok(Self::Header(header))
            }
            Kind::ReferenceSequence => {
                let fields = split_fields(v)?;
                let reference_sequence =
                    Map::<ReferenceSequence>::try_from(fields).map_err(|_| ParseError::Invalid)?;
                Ok(Self::ReferenceSequence(reference_sequence))
            }
            Kind::ReadGroup => {
                let fields = split_fields(v)?;
                let read_group =
                    Map::<ReadGroup>::try_from(fields).map_err(|_| ParseError::Invalid)?;
                Ok(Self::ReadGroup(read_group))
            }
            Kind::Program => {
                let fields = split_fields(v)?;
                let program = Map::<Program>::try_from(fields).map_err(|_| ParseError::Invalid)?;
                Ok(Self::Program(program))
            }
            Kind::Comment => Ok(Self::Comment(v.into())),
        }
    }
}

fn split_fields(s: &str) -> Result<Vec<(String, String)>, ParseError> {
    s.split(DELIMITER).map(split_field).collect()
}

fn split_field(s: &str) -> Result<(String, String), ParseError> {
    const FIELD_DELIMITER: char = ':';

    match s.split_once(FIELD_DELIMITER) {
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
    fn test_from_str() {
        assert!(matches!("@HD\tVN:1.6".parse(), Ok(Record::Header(_))));

        assert!(matches!(
            "@SQ\tSN:sq0\tLN:8".parse(),
            Ok(Record::ReferenceSequence(_))
        ));

        assert!(matches!("@RG\tID:rg0".parse(), Ok(Record::ReadGroup(_))));
        assert!(matches!("@PG\tID:pg0".parse(), Ok(Record::Program(_))));

        assert_eq!(
            "@CO\tnoodles".parse(),
            Ok(Record::Comment(String::from("noodles")))
        );

        assert_eq!("@CO\t".parse(), Ok(Record::Comment(String::new())));

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

        assert_eq!("@CO".parse::<Record>(), Err(ParseError::Invalid));
    }
}
