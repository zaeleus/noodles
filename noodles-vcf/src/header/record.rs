//! VCF header record and components.

pub mod key;
pub(crate) mod parser;
pub mod value;

use indexmap::IndexMap;

pub use self::{key::Key, value::Value};

use std::{error, fmt, str::FromStr};

use super::{AlternativeAllele, Contig, FileFormat, Filter, Format, Info, Meta, Pedigree, Sample};

pub(crate) const PREFIX: &str = "##";

/// A VCF header record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Record {
    /// An `ALT` record.
    AlternativeAllele(AlternativeAllele),
    /// An `assembly` record.
    Assembly(String),
    /// A `contig` record.
    Contig(Contig),
    /// A `fileformat` record.
    FileFormat(FileFormat),
    /// A `FILTER` record.
    Filter(Filter),
    /// A `FORMAT` record.
    Format(Format),
    /// An `INFO` record.
    Info(Info),
    /// A `META` record.
    Meta(Meta),
    /// A `PEDIGREE` record.
    Pedigree(Pedigree),
    /// A `pedigreeDB` record.
    PedigreeDb(String),
    /// A `SAMPLE` record.
    Sample(Sample),
    /// A nonstadard record.
    Other(Key, Value),
}

/// An error returned when a raw VCF header record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => f.write_str("invalid input"),
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (_, (raw_key, raw_value)) = parser::parse(s).map_err(|_| ParseError::Invalid)?;

        let value = match raw_value {
            parser::Value::String(s) => Value::String(s),
            parser::Value::Struct(raw_fields) => {
                let mut fields = IndexMap::new();

                for (k, v) in raw_fields {
                    if fields.insert(k, v).is_some() {
                        return Err(ParseError::Invalid);
                    }
                }

                let id = match fields.shift_remove("ID") {
                    Some(value) => value,
                    None => return Err(ParseError::Invalid),
                };

                Value::Struct(id, fields)
            }
        };

        match Key::from(raw_key) {
            key::FILE_FORMAT => match value {
                Value::String(s) => {
                    let file_format = s.parse().map_err(|_| ParseError::Invalid)?;
                    Ok(Self::FileFormat(file_format))
                }
                _ => Err(ParseError::Invalid),
            },
            key::INFO => match value {
                Value::Struct(id, fields) => {
                    // FIXME
                    let info = Info::try_from_fields(id, fields, Default::default())
                        .map_err(|_| ParseError::Invalid)?;
                    Ok(Self::Info(info))
                }
                _ => Err(ParseError::Invalid),
            },
            key::FILTER => match value {
                Value::Struct(id, fields) => {
                    let filter =
                        Filter::try_from_fields(id, fields).map_err(|_| ParseError::Invalid)?;
                    Ok(Self::Filter(filter))
                }
                _ => Err(ParseError::Invalid),
            },
            key::FORMAT => match value {
                Value::Struct(id, fields) => {
                    // FIXME
                    let format = Format::try_from_fields(id, fields, Default::default())
                        .map_err(|_| ParseError::Invalid)?;
                    Ok(Self::Format(format))
                }
                _ => Err(ParseError::Invalid),
            },
            key::ALTERNATIVE_ALLELE => match value {
                Value::Struct(id, fields) => {
                    let alternative_allele = AlternativeAllele::try_from_fields(id, fields)
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
                Value::Struct(id, fields) => {
                    let contig =
                        Contig::try_from_fields(id, fields).map_err(|_| ParseError::Invalid)?;
                    Ok(Self::Contig(contig))
                }
                _ => Err(ParseError::Invalid),
            },
            key::META => match value {
                Value::Struct(id, fields) => {
                    let meta =
                        Meta::try_from_fields(id, fields).map_err(|_| ParseError::Invalid)?;
                    Ok(Self::Meta(meta))
                }
                _ => Err(ParseError::Invalid),
            },
            key::SAMPLE => match value {
                Value::Struct(id, fields) => {
                    let sample =
                        Sample::try_from_fields(id, fields).map_err(|_| ParseError::Invalid)?;
                    Ok(Self::Sample(sample))
                }
                _ => Err(ParseError::Invalid),
            },
            key::PEDIGREE => match value {
                Value::Struct(id, fields) => {
                    let pedigree =
                        Pedigree::try_from_fields(id, fields).map_err(|_| ParseError::Invalid)?;
                    Ok(Self::Pedigree(pedigree))
                }
                _ => Err(ParseError::Invalid),
            },
            key::PEDIGREE_DB => match value {
                Value::String(s) => Ok(Self::PedigreeDb(s)),
                _ => Err(ParseError::Invalid),
            },
            key => Ok(Self::Other(key, value)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), value::TryFromFieldsError> {
        let line = "##fileformat=VCFv4.3";

        assert_eq!(line.parse(), Ok(Record::FileFormat(FileFormat::new(4, 3))));

        let line =
            r#"##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">"#;

        assert!(matches!(line.parse(), Ok(Record::Info(_))));

        assert!("".parse::<Record>().is_err());

        Ok(())
    }
}
