//! VCF header record and components.

pub mod key;
pub(crate) mod parser;
pub mod value;

pub use self::{key::Key, value::Value};

use std::{error, fmt, str::FromStr};

use self::value::{
    map::{self, AlternativeAllele, Contig, Filter, Format, Info, Meta, Other},
    Map,
};
use super::{file_format, FileFormat};

pub(crate) const PREFIX: &str = "##";

/// A VCF header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Record {
    /// An `ALT` record.
    AlternativeAllele(
        crate::record::alternate_bases::allele::Symbol,
        Map<AlternativeAllele>,
    ),
    /// An `assembly` record.
    Assembly(String),
    /// A `contig` record.
    Contig(map::contig::Name, Map<Contig>),
    /// A `fileformat` record.
    FileFormat(FileFormat),
    /// A `FILTER` record.
    Filter(String, Map<Filter>),
    /// A `FORMAT` record.
    Format(crate::record::genotypes::keys::Key, Map<Format>),
    /// An `INFO` record.
    Info(crate::record::info::field::Key, Map<Info>),
    /// A `META` record.
    Meta(String, Map<Meta>),
    /// A `pedigreeDB` record.
    PedigreeDb(String),
    /// A nonstadard record.
    Other(key::Other, Value),
}

/// An error returned when a raw VCF header record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
    /// The file format record is invalid.
    InvalidFileFormat(file_format::ParseError),
    /// An INFO record is invalid.
    InvalidInfo(
        Option<crate::record::info::field::Key>,
        map::info::ParseError,
    ),
    /// A FILTER record is invalid.
    InvalidFilter(map::filter::ParseError),
    /// A FORMAT record is invalid.
    InvalidFormat(
        Option<crate::record::genotypes::keys::Key>,
        map::format::ParseError,
    ),
    /// An ALT record is invalid.
    InvalidAlternativeAllele(map::alternative_allele::ParseError),
    /// A contig record is invalid.
    InvalidContig(map::contig::ParseError),
    /// A META record is invalid.
    InvalidMeta(map::TryFromFieldsError),
    /// A nonstandard record is invalid.
    InvalidOther(key::Other, map::TryFromFieldsError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Invalid => None,
            Self::InvalidFileFormat(e) => Some(e),
            Self::InvalidInfo(_, e) => Some(e),
            Self::InvalidFormat(_, e) => Some(e),
            Self::InvalidFilter(e) => Some(e),
            Self::InvalidContig(e) => Some(e),
            Self::InvalidAlternativeAllele(e) => Some(e),
            Self::InvalidMeta(e) | Self::InvalidOther(_, e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidFileFormat(_) => write!(f, "invalid {}", key::FILE_FORMAT),
            Self::InvalidInfo(id, _) => {
                write!(f, "invalid {} record", key::INFO)?;

                if let Some(id) = id {
                    write!(f, ": ID={id}")?;
                }

                Ok(())
            }
            Self::InvalidFilter(_) => write!(f, "invalid {} record", key::FILTER),
            Self::InvalidFormat(id, _) => {
                write!(f, "invalid {} record", key::FORMAT)?;

                if let Some(id) = id {
                    write!(f, ": ID={id}")?;
                }

                Ok(())
            }
            Self::InvalidAlternativeAllele(_) => {
                write!(f, "invalid {} record", key::ALTERNATIVE_ALLELE)
            }
            Self::InvalidContig(_) => write!(f, "invalid {} record", key::CONTIG),
            Self::InvalidMeta(_) => write!(f, "invalid {}", key::META),
            Self::InvalidOther(key, _) => write!(f, "invalid {key} record"),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::try_from((FileFormat::default(), s))
    }
}

impl TryFrom<(FileFormat, &str)> for Record {
    type Error = ParseError;

    fn try_from((file_format, s): (FileFormat, &str)) -> Result<Self, Self::Error> {
        const ID: &str = "ID";

        let (_, (raw_key, value)) = parser::parse(s).map_err(|_| ParseError::Invalid)?;

        match Key::from(raw_key) {
            key::FILE_FORMAT => match value {
                parser::Value::String(s) => {
                    let file_format = s.parse().map_err(ParseError::InvalidFileFormat)?;
                    Ok(Self::FileFormat(file_format))
                }
                _ => Err(ParseError::Invalid),
            },
            key::INFO => match value {
                parser::Value::Struct(mut fields) => {
                    let id: crate::record::info::field::Key = remove_field(&mut fields, ID)
                        .ok_or(ParseError::InvalidInfo(
                            None,
                            map::info::ParseError::MissingField(map::info::tag::ID),
                        ))
                        .and_then(|s| {
                            s.parse().map_err(|e| {
                                ParseError::InvalidInfo(None, map::info::ParseError::InvalidId(e))
                            })
                        })?;

                    let info = Map::<Info>::try_from((file_format, fields))
                        .map_err(|e| ParseError::InvalidInfo(Some(id.clone()), e))?;

                    validate_info_definition(file_format, &id, info.number(), info.ty())?;

                    Ok(Self::Info(id, info))
                }
                _ => Err(ParseError::Invalid),
            },
            key::FILTER => match value {
                parser::Value::Struct(mut fields) => {
                    let id = remove_field(&mut fields, ID).ok_or(ParseError::InvalidFilter(
                        map::filter::ParseError::MissingField(map::filter::tag::ID),
                    ))?;

                    let filter =
                        Map::<Filter>::try_from(fields).map_err(ParseError::InvalidFilter)?;

                    Ok(Self::Filter(id, filter))
                }
                _ => Err(ParseError::Invalid),
            },
            key::FORMAT => match value {
                parser::Value::Struct(mut fields) => {
                    let id: crate::record::genotypes::keys::Key = remove_field(&mut fields, ID)
                        .ok_or(ParseError::InvalidFormat(
                            None,
                            map::format::ParseError::MissingField(map::format::tag::ID),
                        ))
                        .and_then(|id| {
                            id.parse().map_err(|e| {
                                ParseError::InvalidFormat(
                                    None,
                                    map::format::ParseError::InvalidId(e),
                                )
                            })
                        })?;

                    let format = Map::<Format>::try_from((file_format, fields))
                        .map_err(|e| ParseError::InvalidFormat(Some(id.clone()), e))?;

                    validate_format_definition(file_format, &id, format.number(), format.ty())?;

                    Ok(Self::Format(id, format))
                }
                _ => Err(ParseError::Invalid),
            },
            key::ALTERNATIVE_ALLELE => match value {
                parser::Value::Struct(mut fields) => {
                    let id = remove_field(&mut fields, ID)
                        .ok_or(ParseError::InvalidAlternativeAllele(
                            map::alternative_allele::ParseError::MissingField(
                                map::alternative_allele::tag::ID,
                            ),
                        ))
                        .and_then(|id| {
                            id.parse().map_err(|e| {
                                ParseError::InvalidAlternativeAllele(
                                    map::alternative_allele::ParseError::InvalidId(e),
                                )
                            })
                        })?;

                    let alternative_allele = Map::<AlternativeAllele>::try_from(fields)
                        .map_err(ParseError::InvalidAlternativeAllele)?;

                    Ok(Self::AlternativeAllele(id, alternative_allele))
                }
                _ => Err(ParseError::Invalid),
            },
            key::ASSEMBLY => match value {
                parser::Value::String(s) => Ok(Self::Assembly(s)),
                _ => Err(ParseError::Invalid),
            },
            key::CONTIG => match value {
                parser::Value::Struct(mut fields) => {
                    let id = remove_field(&mut fields, ID)
                        .ok_or(ParseError::InvalidContig(
                            map::contig::ParseError::MissingField(map::contig::tag::ID),
                        ))
                        .and_then(|id| id.parse().map_err(|_| ParseError::Invalid))?;

                    let contig =
                        Map::<Contig>::try_from(fields).map_err(ParseError::InvalidContig)?;

                    Ok(Self::Contig(id, contig))
                }
                _ => Err(ParseError::Invalid),
            },
            key::META => match value {
                parser::Value::Struct(mut fields) => {
                    let id = remove_field(&mut fields, ID).ok_or(ParseError::InvalidMeta(
                        map::TryFromFieldsError::MissingField(ID),
                    ))?;

                    let meta = Map::<Meta>::try_from(fields).map_err(ParseError::InvalidMeta)?;

                    Ok(Self::Meta(id, meta))
                }
                _ => Err(ParseError::Invalid),
            },
            key::PEDIGREE_DB => match value {
                parser::Value::String(s) => Ok(Self::PedigreeDb(s)),
                _ => Err(ParseError::Invalid),
            },
            Key::Other(k) => {
                let v = match value {
                    parser::Value::String(s) => Value::from(s),
                    parser::Value::Struct(mut fields) => {
                        let id = remove_field(&mut fields, ID).ok_or_else(|| {
                            ParseError::InvalidOther(
                                k.clone(),
                                map::TryFromFieldsError::MissingField(ID),
                            )
                        })?;

                        let map = Map::<Other>::try_from(fields)
                            .map_err(|e| ParseError::InvalidOther(k.clone(), e))?;

                        Value::from((id, map))
                    }
                };

                Ok(Self::Other(k, v))
            }
        }
    }
}

fn remove_field(fields: &mut Vec<(String, String)>, key: &str) -> Option<String> {
    let i = fields.iter().position(|(k, _)| k == key)?;
    let (_, value) = fields.remove(i);
    Some(value)
}

fn validate_format_definition(
    file_format: FileFormat,
    id: &crate::record::genotypes::keys::Key,
    actual_number: super::Number,
    actual_type: super::record::value::map::format::Type,
) -> Result<(), ParseError> {
    use crate::header::record::value::map::format::definition::definition;

    if let Some((expected_number, expected_type, _)) = definition(file_format, id) {
        if actual_number != expected_number {
            return Err(ParseError::InvalidFormat(
                Some(id.clone()),
                map::format::ParseError::NumberMismatch {
                    actual: actual_number,
                    expected: expected_number,
                },
            ));
        }

        if actual_type != expected_type {
            return Err(ParseError::InvalidFormat(
                Some(id.clone()),
                map::format::ParseError::TypeMismatch {
                    actual: actual_type,
                    expected: expected_type,
                },
            ));
        }
    }

    Ok(())
}

fn validate_info_definition(
    file_format: FileFormat,
    id: &crate::record::info::field::Key,
    actual_number: super::Number,
    actual_type: super::record::value::map::info::Type,
) -> Result<(), ParseError> {
    use super::record::value::map::info::definition::definition;

    if let Some((expected_number, expected_type, _)) = definition(file_format, id) {
        if actual_number != expected_number {
            return Err(ParseError::InvalidInfo(
                Some(id.clone()),
                map::info::ParseError::NumberMismatch {
                    actual: actual_number,
                    expected: expected_number,
                },
            ));
        }

        if actual_type != expected_type {
            return Err(ParseError::InvalidInfo(
                Some(id.clone()),
                map::info::ParseError::TypeMismatch {
                    actual: actual_type,
                    expected: expected_type,
                },
            ));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let line = "##fileformat=VCFv4.3";
        assert_eq!(line.parse(), Ok(Record::FileFormat(FileFormat::new(4, 3))));

        let line =
            r#"##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">"#;
        assert!(matches!(line.parse(), Ok(Record::Info(..))));

        assert!("".parse::<Record>().is_err());

        Ok(())
    }
}
