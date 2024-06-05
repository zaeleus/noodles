pub mod alternative_allele;
pub mod contig;
mod field;
pub mod filter;
pub mod format;
pub mod info;
pub mod other;

use std::{error, fmt};

use self::field::split_field;
pub use self::{
    alternative_allele::parse_alternative_allele, contig::parse_contig, filter::parse_filter,
    format::parse_format, info::parse_info, other::parse_other,
};
use crate::header::FileFormat;

const PREFIX: u8 = b'<';
const VCF_4_3: FileFormat = FileFormat::new(4, 3);

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

pub fn is_map(src: &[u8], file_format: FileFormat) -> bool {
    const ID: &[u8] = b"ID=";

    fn contains(buf: &[u8], query: &[u8]) -> bool {
        buf.windows(query.len()).any(|window| window == query)
    }

    let has_prefix = src.first().is_some_and(|&b| b == PREFIX);

    if file_format < VCF_4_3 {
        has_prefix && contains(src, ID)
    } else {
        has_prefix
    }
}

pub fn consume_prefix(src: &mut &[u8]) -> Result<(), ParseError> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_map() {
        const VCF_4_2: FileFormat = FileFormat::new(4, 2);

        assert!(is_map(b"<ID=noodles>", VCF_4_2));
        assert!(!is_map(b"<noodles>", VCF_4_2));
        assert!(!is_map(b"noodles", VCF_4_2));
        assert!(!is_map(b"", VCF_4_2));

        assert!(is_map(b"<ID=noodles>", VCF_4_3));
        assert!(is_map(b"<noodles>", VCF_4_3));
        assert!(!is_map(b"noodles", VCF_4_3));
        assert!(!is_map(b"", VCF_4_3));
    }
}
