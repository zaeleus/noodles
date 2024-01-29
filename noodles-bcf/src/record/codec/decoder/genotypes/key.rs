use std::{error, fmt};

use noodles_vcf as vcf;

use crate::{
    header::string_maps::StringStringMap,
    record::codec::decoder::string_map::{self, read_string_map_entry},
};

pub(super) fn read_key<'h>(
    src: &mut &[u8],
    formats: &'h vcf::header::Formats,
    string_map: &StringStringMap,
) -> Result<&'h str, DecodeError> {
    read_string_map_entry(src, string_map)
        .map_err(DecodeError::InvalidStringMap)
        .and_then(|raw_key| {
            formats
                .get_key_value(raw_key)
                .map(|(k, _)| k.as_ref())
                .ok_or(DecodeError::MissingKey)
        })
}

#[derive(Debug, Eq, PartialEq)]
pub enum DecodeError {
    InvalidStringMap(string_map::DecodeError),
    MissingKey,
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidStringMap(e) => Some(e),
            Self::MissingKey => None,
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidStringMap(_) => write!(f, "invalid string map"),
            Self::MissingKey => write!(f, "missing key"),
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

        let formats = [(String::from(key::GENOTYPE), Map::from(key::GENOTYPE))]
            .into_iter()
            .collect();

        let mut string_map = StringStringMap::default();
        string_map.insert("PASS".into());
        string_map.insert("GT".into());

        assert_eq!(read_key(&mut src, &formats, &string_map), Ok(key::GENOTYPE),);
    }
}
