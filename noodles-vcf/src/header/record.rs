//! VCF header record and components.

pub mod key;
pub(crate) mod parser;
pub mod value;

pub use self::key::Key;

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
    AlternativeAllele(Map<AlternativeAllele>),
    /// An `assembly` record.
    Assembly(String),
    /// A `contig` record.
    Contig(Map<Contig>),
    /// A `fileformat` record.
    FileFormat(FileFormat),
    /// A `FILTER` record.
    Filter(Map<Filter>),
    /// A `FORMAT` record.
    Format(Map<Format>),
    /// An `INFO` record.
    Info(Map<Info>),
    /// A `META` record.
    Meta(Map<Meta>),
    /// A `pedigreeDB` record.
    PedigreeDb(String),
    /// A nonstadard record.
    Other(Key, value::Other),
}

/// An error returned when a raw VCF header record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
    /// The file format record is invalid.
    InvalidFileFormat(file_format::ParseError),
    /// An INFO record is invalid.
    InvalidInfo(map::TryFromFieldsError),
    /// A FILTER record is invalid.
    InvalidFilter(map::TryFromFieldsError),
    /// A FORMAT record is invalid.
    InvalidFormat(map::TryFromFieldsError),
    /// An ALT record is invalid.
    InvalidAlternativeAllele(map::TryFromFieldsError),
    /// A contig record is invalid.
    InvalidContig(map::TryFromFieldsError),
    /// A META record is invalid.
    InvalidMeta(map::TryFromFieldsError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => f.write_str("invalid input"),
            Self::InvalidFileFormat(e) => write!(f, "invalid file format: {}", e),
            Self::InvalidInfo(e) => write!(f, "invalid INFO: {}", e),
            Self::InvalidFilter(e) => write!(f, "invalid FILTER: {}", e),
            Self::InvalidFormat(e) => write!(f, "invalid FORMAT: {}", e),
            Self::InvalidAlternativeAllele(e) => write!(f, "invalid ALT: {}", e),
            Self::InvalidContig(e) => write!(f, "invalid contig: {}", e),
            Self::InvalidMeta(e) => write!(f, "invalid META: {}", e),
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
        use self::parser::Value;

        let (_, (raw_key, value)) = parser::parse(s).map_err(|_| ParseError::Invalid)?;

        match Key::from(raw_key) {
            key::FILE_FORMAT => match value {
                Value::String(s) => {
                    let file_format = s.parse().map_err(ParseError::InvalidFileFormat)?;
                    Ok(Self::FileFormat(file_format))
                }
                _ => Err(ParseError::Invalid),
            },
            key::INFO => match value {
                Value::Struct(fields) => {
                    let info = Map::<Info>::try_from((file_format, fields))
                        .map_err(ParseError::InvalidInfo)?;
                    Ok(Self::Info(info))
                }
                _ => Err(ParseError::Invalid),
            },
            key::FILTER => match value {
                Value::Struct(fields) => {
                    let filter =
                        Map::<Filter>::try_from(fields).map_err(|_| ParseError::Invalid)?;
                    Ok(Self::Filter(filter))
                }
                _ => Err(ParseError::Invalid),
            },
            key::FORMAT => match value {
                Value::Struct(fields) => {
                    let format = Map::<Format>::try_from((file_format, fields))
                        .map_err(|_| ParseError::Invalid)?;
                    Ok(Self::Format(format))
                }
                _ => Err(ParseError::Invalid),
            },
            key::ALTERNATIVE_ALLELE => match value {
                Value::Struct(fields) => {
                    let alternative_allele = Map::<AlternativeAllele>::try_from(fields)
                        .map_err(|_| ParseError::Invalid)?;
                    Ok(Self::AlternativeAllele(alternative_allele))
                }
                _ => Err(ParseError::Invalid),
            },
            key::ASSEMBLY => match value {
                Value::String(s) => Ok(Self::Assembly(s)),
                _ => Err(ParseError::Invalid),
            },
            key::CONTIG => match value {
                Value::Struct(fields) => {
                    let contig =
                        Map::<Contig>::try_from(fields).map_err(|_| ParseError::Invalid)?;
                    Ok(Self::Contig(contig))
                }
                _ => Err(ParseError::Invalid),
            },
            key::META => match value {
                Value::Struct(fields) => {
                    let meta = Map::<Meta>::try_from(fields).map_err(|_| ParseError::Invalid)?;
                    Ok(Self::Meta(meta))
                }
                _ => Err(ParseError::Invalid),
            },
            key::PEDIGREE_DB => match value {
                Value::String(s) => Ok(Self::PedigreeDb(s)),
                _ => Err(ParseError::Invalid),
            },
            k => {
                let v = match value {
                    Value::String(s) => value::Other::from(s),
                    Value::Struct(fields) => {
                        let map =
                            Map::<Other>::try_from(fields).map_err(|_| ParseError::Invalid)?;
                        value::Other::from(map)
                    }
                };

                Ok(Self::Other(k, v))
            }
        }
    }
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
        assert!(matches!(line.parse(), Ok(Record::Info(_))));

        assert!("".parse::<Record>().is_err());

        Ok(())
    }
}
