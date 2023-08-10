mod map;
mod string;

use std::{error, fmt};

use self::string::parse_string;
use crate::header::{
    record::{key, Key},
    Record,
};

/// An error returned when a VCF header record value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidFileFormat,
    InvalidInfo(map::info::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidInfo(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidFileFormat => write!(f, "invalid fileformat"),
            Self::InvalidInfo(_) => write!(f, "invalid info"),
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
        key::FILTER => todo!(),
        key::FORMAT => todo!(),
        key::ALTERNATIVE_ALLELE => todo!(),
        key::ASSEMBLY => todo!(),
        key::CONTIG => todo!(),
        key::META => todo!(),
        key::PEDIGREE_DB => todo!(),
        Key::Other(_) => todo!(),
    }
}
