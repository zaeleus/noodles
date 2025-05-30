use std::{error, fmt};

use noodles_vcf as vcf;

use crate::record::codec::decoder::string_map::{self, read_string_map_entry};

pub(super) fn read_key<'h>(
    src: &mut &[u8],
    header: &'h vcf::Header,
) -> Result<&'h str, DecodeError> {
    read_string_map_entry(src, header.string_maps().strings())
        .map_err(DecodeError::InvalidStringMap)
        .and_then(|raw_key| {
            header
                .formats()
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
    fn test_read_key() -> Result<(), vcf::header::ParseError> {
        use vcf::{
            header::{StringMaps, record::value::Map},
            variant::record::samples::keys::key,
        };

        // Some(Type::Int8(Some(Int8::Value(1))))
        let mut src = &[0x11, 0x01][..];

        let mut header = vcf::Header::builder()
            .add_format(key::GENOTYPE, Map::from(key::GENOTYPE))
            .build();

        *header.string_maps_mut() = StringMaps::try_from(&header)?;

        assert_eq!(read_key(&mut src, &header), Ok(key::GENOTYPE));

        Ok(())
    }
}
