mod key;
mod values;

use std::{error, fmt};

use noodles_vcf::{
    self as vcf,
    record::{genotypes::Keys, Genotypes},
};

use self::{
    key::read_key,
    values::{read_genotype_values, read_values},
};
use crate::header::string_maps::StringStringMap;

pub fn read_genotypes(
    src: &mut &[u8],
    formats: &vcf::header::Formats,
    string_map: &StringStringMap,
    sample_count: usize,
    format_count: usize,
) -> Result<Genotypes, DecodeError> {
    use vcf::record::genotypes::keys::key;

    let mut keys = Vec::with_capacity(format_count);
    let mut samples = vec![Vec::new(); sample_count];

    for _ in 0..format_count {
        let key = read_key(src, formats, string_map).map_err(DecodeError::InvalidKey)?;

        let values = if key == &key::GENOTYPE {
            read_genotype_values(src, sample_count).map_err(DecodeError::InvalidValues)?
        } else {
            read_values(src, sample_count).map_err(DecodeError::InvalidValues)?
        };

        keys.push(key.clone());

        for (sample, value) in samples.iter_mut().zip(values) {
            sample.push(value);
        }
    }

    let keys = Keys::try_from(keys).map_err(DecodeError::InvalidKeys)?;

    Ok(Genotypes::new(keys, samples))
}

#[allow(clippy::enum_variant_names)]
#[derive(Debug, Eq, PartialEq)]
pub enum DecodeError {
    InvalidKey(key::DecodeError),
    InvalidValues(values::DecodeError),
    InvalidKeys(vcf::record::genotypes::keys::TryFromKeyVectorError),
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidKey(e) => Some(e),
            Self::InvalidValues(e) => Some(e),
            Self::InvalidKeys(e) => Some(e),
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidKey(_) => write!(f, "invalid key"),
            Self::InvalidValues(_) => write!(f, "invalid values"),
            Self::InvalidKeys(_) => write!(f, "invalid keys"),
        }
    }
}
