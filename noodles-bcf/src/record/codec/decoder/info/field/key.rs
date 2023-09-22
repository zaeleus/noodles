use std::{error, fmt};

use crate::{
    header::string_maps::StringStringMap, record::codec::decoder::string_map::read_string_map_index,
};

pub(super) fn read_key<'m>(
    src: &mut &[u8],
    string_map: &'m StringStringMap,
) -> Result<&'m str, DecodeError> {
    read_string_map_index(src)
        .map_err(DecodeError::InvalidStringMapIndex)
        .and_then(|j| {
            string_map
                .get_index(j)
                .ok_or(DecodeError::MissingStringMapEntry)
        })
}

#[derive(Debug, Eq, PartialEq)]
pub enum DecodeError {
    InvalidStringMapIndex(crate::record::codec::decoder::string_map::DecodeError),
    MissingStringMapEntry,
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
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_key() {
        // Some(Type::Int8(Some(Int8::Value(1))))
        let mut src = &[0x11, 0x01][..];

        let mut string_map = StringStringMap::default();
        string_map.insert("PASS".into());
        string_map.insert("NS".into());

        assert_eq!(read_key(&mut src, &string_map), Ok("NS"));
    }
}
