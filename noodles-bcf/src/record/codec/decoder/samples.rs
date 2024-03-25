mod key;
mod values;

use std::{error, fmt};

use noodles_vcf::{self as vcf, variant::record_buf::Samples};

use self::{
    key::read_key,
    values::{read_genotype_values, read_values},
};

pub fn read_samples(
    src: &mut &[u8],
    header: &vcf::Header,
    sample_count: usize,
    format_count: usize,
) -> Result<Samples, DecodeError> {
    use vcf::variant::record::samples::keys::key;

    let mut keys = Vec::with_capacity(format_count);
    let mut samples = vec![Vec::new(); sample_count];

    for _ in 0..format_count {
        let key = read_key(src, header).map_err(DecodeError::InvalidKey)?;

        let values = if key == key::GENOTYPE {
            read_genotype_values(src, sample_count).map_err(DecodeError::InvalidValues)?
        } else {
            read_values(src, sample_count).map_err(DecodeError::InvalidValues)?
        };

        keys.push(key.into());

        for (sample, value) in samples.iter_mut().zip(values) {
            sample.push(value);
        }
    }

    Ok(Samples::new(keys.into_iter().collect(), samples))
}

#[allow(clippy::enum_variant_names)]
#[derive(Debug, Eq, PartialEq)]
pub enum DecodeError {
    InvalidKey(key::DecodeError),
    InvalidValues(values::DecodeError),
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidKey(e) => Some(e),
            Self::InvalidValues(e) => Some(e),
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidKey(_) => write!(f, "invalid key"),
            Self::InvalidValues(_) => write!(f, "invalid values"),
        }
    }
}
