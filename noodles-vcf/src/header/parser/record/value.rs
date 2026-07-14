mod map;
mod string;

use std::{error, fmt};

use self::string::parse_string;
use crate::header::{
    FileFormat, Record,
    record::{Key, Value, key},
};

/// An error returned when a VCF header record value fails to parse.
#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidFileFormat(string::file_format::ParseError),
    InvalidInfo(map::info::ParseError),
    InvalidFilter(map::filter::ParseError),
    InvalidFormat(map::format::ParseError),
    InvalidAlternativeAllele(map::alternative_allele::ParseError),
    InvalidContig(map::contig::ParseError),
    InvalidOtherString(key::Other, string::ParseError),
    InvalidOtherMap(key::Other, map::other::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidFileFormat(e) => Some(e),
            Self::InvalidInfo(e) => Some(e),
            Self::InvalidFilter(e) => Some(e),
            Self::InvalidFormat(e) => Some(e),
            Self::InvalidAlternativeAllele(e) => Some(e),
            Self::InvalidContig(e) => Some(e),
            Self::InvalidOtherString(_, e) => Some(e),
            Self::InvalidOtherMap(_, e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidFileFormat(_) => write!(f, "invalid fileformat"),
            Self::InvalidInfo(e) => {
                write!(f, "invalid {}", key::INFO)?;

                if let Some(id) = e.id() {
                    write!(f, ": ID={id}")?;
                }

                Ok(())
            }
            Self::InvalidFilter(e) => {
                write!(f, "invalid {}", key::FILTER)?;

                if let Some(id) = e.id() {
                    write!(f, ": ID={id}")?;
                }

                Ok(())
            }
            Self::InvalidFormat(e) => {
                write!(f, "invalid {}", key::FORMAT)?;

                if let Some(id) = e.id() {
                    write!(f, ": ID={id}")?;
                }

                Ok(())
            }
            Self::InvalidAlternativeAllele(e) => {
                write!(
                    f,
                    "invalid alternative allele ({})",
                    key::ALTERNATIVE_ALLELE
                )?;

                if let Some(id) = e.id() {
                    write!(f, ": ID={id}")?;
                }

                Ok(())
            }
            Self::InvalidContig(e) => {
                write!(f, "invalid {}", key::CONTIG)?;

                if let Some(id) = e.id() {
                    write!(f, ": ID={id}")?;
                }

                Ok(())
            }
            Self::InvalidOtherString(key, _) => write!(f, "invalid other string: {key}"),
            Self::InvalidOtherMap(key, e) => {
                write!(f, "invalid other map: {key}")?;

                if let Some(id) = e.id() {
                    write!(f, ": ID={id}")?;
                }

                Ok(())
            }
        }
    }
}

pub(super) fn parse_value(
    src: &mut &[u8],
    file_format: FileFormat,
    key: Key,
) -> Result<Record, ParseError> {
    const META: &str = "META";
    const PEDIGREE: &str = "PEDIGREE";

    match key {
        key::FILE_FORMAT => string::parse_file_format(src)
            .map(Record::FileFormat)
            .map_err(ParseError::InvalidFileFormat),
        key::INFO => map::parse_info(src, file_format)
            .map(|(id, map)| Record::Info(id, map))
            .map_err(ParseError::InvalidInfo),
        key::FILTER => map::parse_filter(src)
            .map(|(id, map)| Record::Filter(id, map))
            .map_err(ParseError::InvalidFilter),
        key::FORMAT => map::parse_format(src, file_format)
            .map(|(id, map)| Record::Format(id, map))
            .map_err(ParseError::InvalidFormat),
        key::ALTERNATIVE_ALLELE => map::parse_alternative_allele(src)
            .map(|(id, map)| Record::AlternativeAllele(id, map))
            .map_err(ParseError::InvalidAlternativeAllele),
        key::CONTIG => map::parse_contig(src)
            .map(|(id, map)| Record::Contig(id, map))
            .map_err(ParseError::InvalidContig),
        Key::Other(k) => {
            let v = if k.as_ref() == META {
                map::other::parse_meta(src, file_format)
                    .map(Value::from)
                    .map_err(|e| ParseError::InvalidOtherMap(k.clone(), e))?
            } else if k.as_ref() == PEDIGREE {
                map::other::parse_pedigree(src, file_format)
                    .map(Value::from)
                    .map_err(|e| ParseError::InvalidOtherMap(k.clone(), e))?
            } else if map::is_map(src, file_format) {
                map::parse_other(src)
                    .map(Value::from)
                    .map_err(|e| ParseError::InvalidOtherMap(k.clone(), e))?
            } else {
                parse_string(src)
                    .map(Value::from)
                    .map_err(|e| ParseError::InvalidOtherString(k.clone(), e))?
            };

            Ok(Record::Other(k, v))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_value_with_redefined_reserved_format_key() {
        use crate::header::record::value::map::format::{Number, Type};

        // A header may loosen a reserved field's `Number` (e.g. `DP` from the reserved
        // `Number=1` to `Number=.`). Such a definition is valid VCF and must be honored
        // rather than rejected against the version-reserved default.
        let mut src = &b"<ID=DP,Number=.,Type=Integer,Description=\"d\">"[..];
        let record = parse_value(&mut src, FileFormat::new(4, 4), key::FORMAT).unwrap();

        match record {
            Record::Format(id, map) => {
                assert_eq!(id, "DP");
                assert_eq!(map.number(), Number::Unknown);
                assert_eq!(map.ty(), Type::Integer);
            }
            _ => panic!("expected a FORMAT record"),
        }

        // `AD` is reserved as `Number=R`; a header declaring `Number=.` must be honored.
        let mut src = &b"<ID=AD,Number=.,Type=Integer,Description=\"a\">"[..];
        let record = parse_value(&mut src, FileFormat::new(4, 4), key::FORMAT).unwrap();

        match record {
            Record::Format(id, map) => {
                assert_eq!(id, "AD");
                assert_eq!(map.number(), Number::Unknown);
                assert_eq!(map.ty(), Type::Integer);
            }
            _ => panic!("expected a FORMAT record"),
        }
    }

    #[test]
    fn test_parse_value_with_redefined_reserved_info_key() {
        use crate::header::record::value::map::info::{Number, Type};

        // `DP` is reserved in INFO as `Number=1`; a header declaring `Number=.` must be
        // honored rather than rejected.
        let mut src = &b"<ID=DP,Number=.,Type=Integer,Description=\"d\">"[..];
        let record = parse_value(&mut src, FileFormat::new(4, 4), key::INFO).unwrap();

        match record {
            Record::Info(id, map) => {
                assert_eq!(id, "DP");
                assert_eq!(map.number(), Number::Unknown);
                assert_eq!(map.ty(), Type::Integer);
            }
            _ => panic!("expected an INFO record"),
        }
    }
}
