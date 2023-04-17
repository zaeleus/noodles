//! VCF record genotypes keys.

pub mod key;

pub use self::key::Key;

use std::{
    error, fmt,
    ops::{Deref, DerefMut},
    str::FromStr,
};

use indexmap::IndexSet;

use crate::{header, record::MISSING_FIELD};

const DELIMITER: char = ':';

/// A VCF record genotypes keys, i.e., `FORMAT`.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Keys(IndexSet<Key>);

impl Keys {
    /// Parses raw VCF record genotypes keys.
    pub fn try_from_str(s: &str, formats: &header::Formats) -> Result<Self, ParseError> {
        parse(s, formats)
    }
}

impl Deref for Keys {
    type Target = IndexSet<Key>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Keys {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl fmt::Display for Keys {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            f.write_str(MISSING_FIELD)
        } else {
            for (i, key) in self.iter().enumerate() {
                if i > 0 {
                    write!(f, "{DELIMITER}")?;
                }

                f.write_str(key.as_ref())?;
            }

            Ok(())
        }
    }
}

/// An error returned when raw VCF record genotypes keys fail to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The key is invalid.
    InvalidKey(key::ParseError),
    /// The format is invalid.
    InvalidFormat(TryFromKeyVectorError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Empty => None,
            Self::InvalidKey(e) => Some(e),
            Self::InvalidFormat(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::InvalidKey(_) => f.write_str("invalid key"),
            Self::InvalidFormat(_) => f.write_str("invalid format"),
        }
    }
}

impl FromStr for Keys {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse(s, &header::Formats::default())
    }
}

fn parse(s: &str, formats: &header::Formats) -> Result<Keys, ParseError> {
    if s.is_empty() {
        Err(ParseError::Empty)
    } else {
        s.split(DELIMITER)
            .map(
                |raw_key| match formats.keys().find(|k| k.as_ref() == raw_key) {
                    Some(k) => Ok(k.clone()),
                    None => raw_key.parse(),
                },
            )
            .collect::<Result<Vec<_>, _>>()
            .map_err(ParseError::InvalidKey)
            .and_then(|keys| Keys::try_from(keys).map_err(ParseError::InvalidFormat))
    }
}

/// An error returned when a vector of record genotypes keys fails to convert to a format.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromKeyVectorError {
    /// The genotype key (`GT`) position is invalid.
    ///
    /// The genotype key must be first if present. See § 1.6.2 Genotype fields (2020-06-25).
    InvalidGenotypeKeyPosition,
    /// A key is duplicated.
    ///
    /// § 1.6.2 Genotype fields (2021-01-13): "...duplicate keys are not allowed".
    DuplicateKey(Key),
}

impl error::Error for TryFromKeyVectorError {}

impl fmt::Display for TryFromKeyVectorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidGenotypeKeyPosition => f.write_str("invalid genotype key position"),
            Self::DuplicateKey(key) => write!(f, "duplicate key: {key}"),
        }
    }
}

impl TryFrom<Vec<Key>> for Keys {
    type Error = TryFromKeyVectorError;

    fn try_from(keys: Vec<Key>) -> Result<Self, Self::Error> {
        if keys.is_empty() {
            return Ok(Keys::default());
        } else if let Some(i) = keys.iter().position(|k| k == &key::GENOTYPE) {
            if i != 0 {
                return Err(TryFromKeyVectorError::InvalidGenotypeKeyPosition);
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
    fn test_fmt() {
        let keys = Keys::default();
        assert_eq!(keys.to_string(), ".");

        let keys = Keys([key::GENOTYPE].into_iter().collect());
        assert_eq!(keys.to_string(), "GT");

        let keys = Keys(
            [
                key::GENOTYPE,
                key::CONDITIONAL_GENOTYPE_QUALITY,
                key::READ_DEPTH,
                key::HAPLOTYPE_QUALITY,
            ]
            .into_iter()
            .collect(),
        );
        assert_eq!(keys.to_string(), "GT:GQ:DP:HQ");
    }

    #[test]
    fn test_from_str() {
        assert_eq!(
            "GT".parse(),
            Ok(Keys([key::GENOTYPE].into_iter().collect()))
        );
        assert_eq!(
            "GT:GQ".parse(),
            Ok(Keys(
                [key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY]
                    .into_iter()
                    .collect()
            ))
        );

        assert_eq!("".parse::<Keys>(), Err(ParseError::Empty));
        assert!(matches!(
            "GQ:GT".parse::<Keys>(),
            Err(ParseError::InvalidFormat(_))
        ));
    }

    #[test]
    fn test_try_from_vec_key_for_format() {
        assert_eq!(Keys::try_from(Vec::new()), Ok(Keys::default()));

        assert_eq!(
            Keys::try_from(vec![key::GENOTYPE]),
            Ok(Keys([key::GENOTYPE].into_iter().collect()))
        );

        assert_eq!(
            Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY]),
            Ok(Keys(
                [key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY]
                    .into_iter()
                    .collect()
            ))
        );

        assert_eq!(
            Keys::try_from(vec![key::CONDITIONAL_GENOTYPE_QUALITY]),
            Ok(Keys(
                [key::CONDITIONAL_GENOTYPE_QUALITY].into_iter().collect()
            ))
        );

        assert_eq!(
            Keys::try_from(vec![key::CONDITIONAL_GENOTYPE_QUALITY, key::GENOTYPE]),
            Err(TryFromKeyVectorError::InvalidGenotypeKeyPosition)
        );

        assert_eq!(
            Keys::try_from(vec![key::GENOTYPE, key::GENOTYPE]),
            Err(TryFromKeyVectorError::DuplicateKey(key::GENOTYPE))
        );
    }
}
