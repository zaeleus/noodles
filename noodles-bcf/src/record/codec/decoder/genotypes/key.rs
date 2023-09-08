use std::{error, fmt};

use noodles_vcf::{self as vcf, record::genotypes::keys::Key};

use crate::{
    header::string_maps::StringStringMap, record::codec::decoder::string_map::read_string_map_index,
};

pub(super) fn read_key<'h>(
    src: &mut &[u8],
    formats: &'h vcf::header::Formats,
    string_map: &StringStringMap,
) -> Result<&'h Key, DecodeError> {
    read_string_map_index(src)
        .map_err(DecodeError::InvalidStringMapIndex)
        .and_then(|j| {
            string_map
                .get_index(j)
                .ok_or(DecodeError::MissingStringMapEntry)
        })
        .and_then(|raw_key| {
            formats
                .get_key_value(raw_key)
                .map(|(k, _)| k)
                .ok_or(DecodeError::MissingFormatMapEntry)
        })
}

#[derive(Debug, Eq, PartialEq)]
pub enum DecodeError {
    InvalidStringMapIndex(crate::record::codec::decoder::string_map::DecodeError),
    MissingStringMapEntry,
    MissingFormatMapEntry,
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidStringMapIndex(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidStringMapIndex(_) => write!(f, "invalid string map index"),
            Self::MissingStringMapEntry => write!(f, "missing string map entry"),
            Self::MissingFormatMapEntry => write!(f, "missing format map entry"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_key() {
        use vcf::{header::record::value::Map, record::genotypes::keys::key};

        // Some(Type::Int8(Some(Int8::Value(1))))
        let mut src = &[0x11, 0x01][..];

        let formats = [(key::GENOTYPE, Map::from(&key::GENOTYPE))]
            .into_iter()
            .collect();

        let mut string_map = StringStringMap::default();
        string_map.insert("PASS".into());
        string_map.insert("GT".into());

        assert_eq!(
            read_key(&mut src, &formats, &string_map),
            Ok(&key::GENOTYPE),
        );
    }
}
