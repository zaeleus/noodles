//! VCF record genotypes keys.

pub mod key;

use std::{
    error, fmt,
    ops::{Deref, DerefMut},
};

use indexmap::IndexSet;

/// A VCF record genotypes keys, i.e., `FORMAT`.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Keys(IndexSet<String>);

impl Deref for Keys {
    type Target = IndexSet<String>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Keys {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

/// An error returned when raw VCF record genotypes keys fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The format is invalid.
    InvalidFormat(TryFromKeyVectorError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Empty => None,
            Self::InvalidFormat(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::InvalidFormat(_) => f.write_str("invalid format"),
        }
    }
}

/// An error returned when a vector of record genotypes keys fails to convert to a format.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromKeyVectorError {
    /// The genotype key (`GT`) position is invalid.
    ///
    /// The genotype key must be first if present. See § 1.6.2 Genotype fields (2020-06-25).
    InvalidGenotypeKeyPosition(usize),
    /// A key is duplicated.
    ///
    /// § 1.6.2 Genotype fields (2021-01-13): "...duplicate keys are not allowed".
    DuplicateKey(String),
}

impl error::Error for TryFromKeyVectorError {}

impl fmt::Display for TryFromKeyVectorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidGenotypeKeyPosition(i) => {
                write!(f, "invalid genotype key position: expected 0, got {i}")
            }
            Self::DuplicateKey(key) => write!(f, "duplicate key: {key}"),
        }
    }
}

impl TryFrom<Vec<String>> for Keys {
    type Error = TryFromKeyVectorError;

    fn try_from(keys: Vec<String>) -> Result<Self, Self::Error> {
        if keys.is_empty() {
            return Ok(Keys::default());
        } else if let Some(i) = keys.iter().position(|k| k == key::GENOTYPE) {
            if i != 0 {
                return Err(TryFromKeyVectorError::InvalidGenotypeKeyPosition(i));
            }
        }

        let mut set = IndexSet::new();

        for key in &keys {
            if !set.insert(key.clone()) {
                return Err(TryFromKeyVectorError::DuplicateKey(key.clone()));
            }
        }

        Ok(Self(set))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_vec_key_for_format() {
        assert_eq!(Keys::try_from(Vec::new()), Ok(Keys::default()));

        assert_eq!(
            Keys::try_from(vec![String::from(key::GENOTYPE)]),
            Ok(Keys(
                [key::GENOTYPE].into_iter().map(String::from).collect()
            ))
        );

        assert_eq!(
            Keys::try_from(vec![
                String::from(key::GENOTYPE),
                String::from(key::CONDITIONAL_GENOTYPE_QUALITY)
            ]),
            Ok(Keys(
                [key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY]
                    .into_iter()
                    .map(String::from)
                    .collect()
            ))
        );

        assert_eq!(
            Keys::try_from(vec![String::from(key::CONDITIONAL_GENOTYPE_QUALITY)]),
            Ok(Keys(
                [key::CONDITIONAL_GENOTYPE_QUALITY]
                    .into_iter()
                    .map(String::from)
                    .collect()
            ))
        );

        assert_eq!(
            Keys::try_from(vec![
                String::from(key::CONDITIONAL_GENOTYPE_QUALITY),
                String::from(key::GENOTYPE)
            ]),
            Err(TryFromKeyVectorError::InvalidGenotypeKeyPosition(1))
        );

        assert_eq!(
            Keys::try_from(vec![
                String::from(key::GENOTYPE),
                String::from(key::GENOTYPE)
            ]),
            Err(TryFromKeyVectorError::DuplicateKey(String::from(
                key::GENOTYPE
            )))
        );
    }
}
