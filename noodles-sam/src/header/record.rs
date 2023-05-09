//! SAM header record.

pub mod kind;
pub mod value;

use std::{error, fmt, str::FromStr};

pub use self::kind::Kind;

use self::value::{
    map::{self, header::Version, Program, ReadGroup, ReferenceSequence},
    Map,
};

const DELIMITER: char = '\t';

/// A SAM header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Record {
    /// A header (`HD`) record.
    Header(Map<map::Header>),
    /// A reference sequence (`SQ`) record.
    ReferenceSequence(map::reference_sequence::Name, Map<ReferenceSequence>),
    /// A read group (`RG`) record.
    ReadGroup(String, Map<ReadGroup>),
    /// A program (`PG`) record.
    Program(String, Map<Program>),
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
    /// A value is invalid.
    InvalidValue,
    /// A header record is invalid.
    InvalidHeader(map::header::ParseError),
    /// A reference sequence record is invalid.
    InvalidReferenceSequence(
        Option<map::reference_sequence::Name>,
        map::reference_sequence::ParseError,
    ),
    /// A read group record is invalid.
    InvalidReadGroup(Option<String>, map::read_group::ParseError),
    /// A program record is invalid.
    InvalidProgram(Option<String>, map::program::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidKind(e) => Some(e),
            Self::InvalidHeader(e) => Some(e),
            Self::InvalidReferenceSequence(_, e) => Some(e),
            Self::InvalidReadGroup(_, e) => Some(e),
            Self::InvalidProgram(_, e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => write!(f, "invalid input"),
            Self::InvalidKind(_) => f.write_str("invalid kind"),
            Self::InvalidField => write!(f, "invalid field"),
            Self::InvalidValue => write!(f, "invalid value"),
            Self::InvalidHeader(_) => f.write_str("invalid header (HD) record"),
            Self::InvalidReferenceSequence(name, _) => {
                write!(f, "invalid reference sequence (SQ) record")?;

                if let Some(name) = name {
                    write!(f, ": {}:{}", map::reference_sequence::tag::NAME, name)?;
                }

                Ok(())
            }
            Self::InvalidReadGroup(id, _) => {
                write!(f, "invalid read group (RG) record")?;

                if let Some(id) = id {
                    write!(f, ": {}:{}", map::read_group::tag::ID, id)?;
                }

                Ok(())
            }
            Self::InvalidProgram(id, _) => {
                write!(f, "invalid program (PG) record")?;

                if let Some(id) = id {
                    write!(f, ": {}:{}", map::program::tag::ID, id)?;
                }

                Ok(())
            }
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        const ID: &str = "ID";
        const SN: &str = "SN";

        let (kind, v) = split_record(s)?;

        match kind {
            Kind::Header => {
                let fields = split_fields(v)?;
                let header = Map::try_from(fields).map_err(ParseError::InvalidHeader)?;
                Ok(Self::Header(header))
            }
            Kind::ReferenceSequence => {
                let mut fields = split_fields(v)?;

                let name: map::reference_sequence::Name = remove_field(&mut fields, SN)
                    .ok_or(ParseError::InvalidReferenceSequence(
                        None,
                        map::reference_sequence::ParseError::MissingField(
                            map::reference_sequence::tag::NAME,
                        ),
                    ))
                    .and_then(|t| {
                        t.parse().map_err(|e| {
                            ParseError::InvalidReferenceSequence(
                                None,
                                map::reference_sequence::ParseError::InvalidName(e),
                            )
                        })
                    })?;

                let reference_sequence = Map::try_from(fields)
                    .map_err(|e| ParseError::InvalidReferenceSequence(Some(name.clone()), e))?;

                Ok(Self::ReferenceSequence(name, reference_sequence))
            }
            Kind::ReadGroup => {
                let mut fields = split_fields(v)?;

                let id = remove_field(&mut fields, ID).ok_or(ParseError::InvalidReadGroup(
                    None,
                    map::read_group::ParseError::MissingField(map::read_group::tag::ID),
                ))?;

                let read_group = Map::try_from(fields)
                    .map_err(|e| ParseError::InvalidReadGroup(Some(id.clone()), e))?;

                Ok(Self::ReadGroup(id, read_group))
            }
            Kind::Program => {
                let mut fields = split_fields(v)?;

                let id = remove_field(&mut fields, ID).ok_or(ParseError::InvalidProgram(
                    None,
                    map::program::ParseError::MissingField(map::program::tag::ID),
                ))?;

                let program = Map::try_from(fields)
                    .map_err(|e| ParseError::InvalidProgram(Some(id.clone()), e))?;

                Ok(Self::Program(id, program))
            }
            Kind::Comment => Ok(Self::Comment(v.into())),
        }
    }
}

pub(super) fn extract_version(s: &str) -> Option<Result<Version, ParseError>> {
    const VN: &str = "VN";

    let (kind, fields) = match split_record(s) {
        Ok((k, v)) => (k, v),
        Err(e) => return Some(Err(e)),
    };

    if kind == Kind::Header {
        for result in fields.split(DELIMITER).map(split_field) {
            let (tag, value) = match result {
                Ok((t, v)) => (t, v),
                Err(e) => return Some(Err(e)),
            };

            if tag == VN {
                return Some(
                    value
                        .parse()
                        .map_err(map::header::ParseError::InvalidVersion)
                        .map_err(ParseError::InvalidHeader),
                );
            }
        }
    }

    None
}

fn split_record(s: &str) -> Result<(Kind, &str), ParseError> {
    let (k, v) = s.split_once(DELIMITER).ok_or(ParseError::Invalid)?;
    let kind = k.parse().map_err(ParseError::InvalidKind)?;
    Ok((kind, v))
}

fn split_fields(s: &str) -> Result<Vec<(String, String)>, ParseError> {
    s.split(DELIMITER)
        .map(split_field)
        .map(|result| result.map(|(key, value)| (key.into(), value.into())))
        .collect()
}

fn split_field(s: &str) -> Result<(&str, &str), ParseError> {
    const FIELD_DELIMITER: char = ':';

    match s.split_once(FIELD_DELIMITER) {
        Some((tag, value)) => {
            if !is_valid_value(value) {
                return Err(ParseError::InvalidValue);
            }

            Ok((tag, value))
        }
        None => Err(ParseError::InvalidField),
    }
}

fn is_valid_value(s: &str) -> bool {
    !s.is_empty() && s.chars().all(|c| matches!(c, ' '..='~'))
}

fn remove_field(fields: &mut Vec<(String, String)>, tag: &str) -> Option<String> {
    let i = fields.iter().position(|(t, _)| t == tag)?;
    let (_, value) = fields.remove(i);
    Some(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert!(matches!("@HD\tVN:1.6".parse(), Ok(Record::Header(_))));

        assert!(matches!(
            "@SQ\tSN:sq0\tLN:8".parse(),
            Ok(Record::ReferenceSequence(..))
        ));

        assert!(matches!(
            "@RG\tID:rg0".parse(),
            Ok(Record::ReadGroup(id, _)) if id == "rg0"
        ));

        assert!(matches!(
            "@PG\tID:pg0".parse(),
            Ok(Record::Program(id, _)) if id == "pg0"
        ));

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

        assert_eq!("@PG\tID:".parse::<Record>(), Err(ParseError::InvalidValue));
        assert_eq!(
            "@PG\tID:🍜".parse::<Record>(),
            Err(ParseError::InvalidValue)
        );

        assert_eq!("@CO".parse::<Record>(), Err(ParseError::Invalid));
    }

    #[test]
    fn test_extract_version() {
        assert_eq!(extract_version("@HD\tVN:1.6"), Some(Ok(Version::new(1, 6))));
        assert_eq!(
            extract_version("@HD\tSO:coordinate\tVN:1.6"),
            Some(Ok(Version::new(1, 6)))
        );
        assert_eq!(extract_version("@SQ\tSN:sq0\tLN:8\tVN:1.6"), None);
        assert_eq!(extract_version("@CO\tVN:1.6"), None);
    }
}
