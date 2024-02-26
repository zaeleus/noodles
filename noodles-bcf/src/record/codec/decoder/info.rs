mod field;

use std::{error, fmt};

use noodles_vcf::{self as vcf, header::string_maps::StringStringMap};

pub(crate) use self::field::read_field;

pub fn read_info(
    src: &mut &[u8],
    infos: &vcf::header::Infos,
    string_string_map: &StringStringMap,
    len: usize,
    info: &mut vcf::variant::record_buf::Info,
) -> Result<(), DecodeError> {
    info.clear();

    for _ in 0..len {
        let (key, value) =
            read_field(src, infos, string_string_map).map_err(DecodeError::InvalidField)?;

        if info.insert(key.clone(), value).is_some() {
            return Err(DecodeError::DuplicateKey(key));
        }
    }

    Ok(())
}

#[derive(Debug, Eq, PartialEq)]
pub enum DecodeError {
    InvalidField(field::DecodeError),
    DuplicateKey(String),
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidField(e) => Some(e),
            Self::DuplicateKey(_) => None,
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField(_) => write!(f, "invalid field"),
            Self::DuplicateKey(key) => write!(f, "duplicate key: {key}"),
        }
    }
}
