mod key;
mod value;

use std::{error, fmt};

use noodles_vcf as vcf;

use self::{key::read_key, value::read_value};
use crate::header::string_maps::StringStringMap;

pub(crate) fn read_field(
    src: &mut &[u8],
    infos: &vcf::header::Infos,
    string_map: &StringStringMap,
) -> Result<
    (
        vcf::record::info::field::Key,
        Option<vcf::record::info::field::Value>,
    ),
    DecodeError,
> {
    let key = read_key(src, infos, string_map).map_err(DecodeError::InvalidKey)?;
    let info = infos.get(key).ok_or(DecodeError::MissingInfoMapEntry)?;
    let value = read_value(src, info).map_err(DecodeError::InvalidValue)?;
    Ok((key.clone(), value))
}

#[derive(Debug, Eq, PartialEq)]
pub enum DecodeError {
    InvalidKey(key::DecodeError),
    MissingInfoMapEntry,
    InvalidValue(value::DecodeError),
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidKey(e) => Some(e),
            Self::MissingInfoMapEntry => None,
            Self::InvalidValue(e) => Some(e),
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidKey(_) => write!(f, "invalid key"),
            Self::MissingInfoMapEntry => write!(f, "missing info map entry"),
            Self::InvalidValue(_) => write!(f, "invalid value"),
        }
    }
}
