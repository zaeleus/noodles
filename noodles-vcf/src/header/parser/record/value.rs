mod map;
mod string;

use std::{error, fmt};

use self::string::parse_string;
use crate::header::{
    record::{key, Key, Value},
    Record,
};

/// An error returned when a VCF header record value fails to parse.
#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidFileFormat,
    InvalidInfo(map::info::ParseError),
    InvalidFilter(map::filter::ParseError),
    InvalidFormat(map::format::ParseError),
    InvalidAlternativeAllele(map::alternative_allele::ParseError),
    InvalidContig(map::contig::ParseError),
    InvalidOther,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidInfo(e) => Some(e),
            Self::InvalidFilter(e) => Some(e),
            Self::InvalidFormat(e) => Some(e),
            Self::InvalidAlternativeAllele(e) => Some(e),
            Self::InvalidContig(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidFileFormat => write!(f, "invalid fileformat"),
            Self::InvalidInfo(_) => write!(f, "invalid info"),
            Self::InvalidFilter(_) => write!(f, "invalid filter"),
            Self::InvalidFormat(_) => write!(f, "invalid format"),
            Self::InvalidAlternativeAllele(_) => write!(f, "invalid alternative allele"),
            Self::InvalidContig(_) => write!(f, "invalid contig"),
            Self::InvalidOther => write!(f, "invalid other"),
        }
    }
}

pub(super) fn parse_value(src: &mut &[u8], key: Key) -> Result<Record, ParseError> {
    match key {
        key::FILE_FORMAT => parse_string(src)
            .map_err(|_| ParseError::InvalidFileFormat)
            .and_then(|s| s.parse().map_err(|_| ParseError::InvalidFileFormat))
            .map(Record::FileFormat),
        key::INFO => map::parse_info(src)
            .map(|(id, map)| Record::Info(id, map))
            .map_err(ParseError::InvalidInfo),
        key::FILTER => map::parse_filter(src)
            .map(|(id, map)| Record::Filter(id, map))
            .map_err(ParseError::InvalidFilter),
        key::FORMAT => map::parse_format(src)
            .map(|(id, map)| Record::Format(id, map))
            .map_err(ParseError::InvalidFormat),
        key::ALTERNATIVE_ALLELE => map::parse_alternative_allele(src)
            .map(|(id, map)| Record::AlternativeAllele(id, map))
            .map_err(ParseError::InvalidAlternativeAllele),
        key::CONTIG => map::parse_contig(src)
            .map(|(id, map)| Record::Contig(id, map))
            .map_err(ParseError::InvalidContig),
        Key::Other(k) => {
            let is_map = src.first().map(|&b| b == b'<').unwrap_or_default();

            let v = if is_map {
                map::parse_other(src)
                    .map(Value::from)
                    .map_err(|_| ParseError::InvalidOther)?
            } else {
                parse_string(src)
                    .map(String::from)
                    .map(Value::from)
                    .map_err(|_| ParseError::InvalidOther)?
            };

            Ok(Record::Other(k, v))
        }
    }
}
