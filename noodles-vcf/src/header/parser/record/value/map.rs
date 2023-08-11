pub mod alternative_allele;
pub mod contig;
mod field;
pub mod filter;
pub mod format;
pub mod info;
pub mod other;

use std::{error, fmt};

pub use self::{
    alternative_allele::parse_alternative_allele, contig::parse_contig, filter::parse_filter,
    format::parse_format, info::parse_info, other::parse_other,
};

/// An error returned when a VCF header record map value fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    InvalidPrefix,
    InvalidSuffix,
    UnexpectedEof,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPrefix => write!(f, "invalid prefix"),
            Self::InvalidSuffix => write!(f, "invalid suffix"),
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
        }
    }
}

pub fn consume_prefix(src: &mut &[u8]) -> Result<(), ParseError> {
    const PREFIX: u8 = b'<';

    if let Some((b, rest)) = src.split_first() {
        if *b == PREFIX {
            *src = rest;
            Ok(())
        } else {
            Err(ParseError::InvalidPrefix)
        }
    } else {
        Err(ParseError::UnexpectedEof)
    }
}

pub fn consume_suffix(src: &mut &[u8]) -> Result<(), ParseError> {
    const SUFFIX: u8 = b'>';

    if let Some((b, rest)) = src.split_first() {
        if *b == SUFFIX {
            *src = rest;
            Ok(())
        } else {
            Err(ParseError::InvalidSuffix)
        }
    } else {
        Err(ParseError::UnexpectedEof)
    }
}
