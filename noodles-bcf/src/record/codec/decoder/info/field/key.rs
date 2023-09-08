use std::{error, fmt};

use noodles_vcf as vcf;

use crate::{
    header::string_maps::StringStringMap, record::codec::decoder::string_map::read_string_map_index,
};

pub(super) fn read_key<'h>(
    src: &mut &[u8],
    infos: &'h vcf::header::Infos,
    string_map: &StringStringMap,
) -> Result<&'h vcf::record::info::field::Key, DecodeError> {
    read_string_map_index(src)
        .map_err(DecodeError::InvalidStringMapIndex)
        .and_then(|j| {
            string_map
                .get_index(j)
                .ok_or(DecodeError::MissingStringMapEntry)
        })
        .and_then(|raw_key| {
            infos
                .get_key_value(raw_key)
                .map(|(k, _)| k)
                .ok_or(DecodeError::MissingInfoMapEntry)
        })
}

#[derive(Debug, Eq, PartialEq)]
pub enum DecodeError {
    InvalidStringMapIndex(crate::record::codec::decoder::string_map::DecodeError),
    MissingStringMapEntry,
    MissingInfoMapEntry,
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
            Self::MissingInfoMapEntry => write!(f, "missing info map entry"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_key() {
        use vcf::{header::record::value::Map, record::info::field::key};

        // Some(Type::Int8(Some(Int8::Value(1))))
        let mut src = &[0x11, 0x01][..];

        let infos = [(
            key::SAMPLES_WITH_DATA_COUNT,
            Map::from(&key::SAMPLES_WITH_DATA_COUNT),
        )]
        .into_iter()
        .collect();

        let mut string_map = StringStringMap::default();
        string_map.insert("PASS".into());
        string_map.insert("NS".into());

        assert_eq!(
            read_key(&mut src, &infos, &string_map),
            Ok(&key::SAMPLES_WITH_DATA_COUNT),
        );
    }
}
