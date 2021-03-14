//! VCF record genotype format.

use std::{collections::HashSet, convert::TryFrom, error, fmt, ops::Deref, str::FromStr};

use super::genotype::field::{key, Key};

const DELIMITER: char = ':';

/// A VCF record genotype format (`FORMAT`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Format(Vec<Key>);

impl Deref for Format {
    type Target = [Key];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for Format {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, key) in self.iter().enumerate() {
            if i > 0 {
                write!(f, "{}", DELIMITER)?
            }

            f.write_str(key.as_ref())?;
        }

        Ok(())
    }
}

/// An error returned when a raw VCF record format fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The key is invalid.
    InvalidKey(key::ParseError),
    /// The format is invalid.
    InvalidFormat(TryFromKeyVectorError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::InvalidKey(e) => write!(f, "{}", e),
            Self::InvalidFormat(e) => write!(f, "{}", e),
        }
    }
}

impl FromStr for Format {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            Err(ParseError::Empty)
        } else {
            s.split(DELIMITER)
                .map(|s| s.parse())
                .collect::<Result<Vec<_>, _>>()
                .map_err(ParseError::InvalidKey)
                .and_then(|keys| Self::try_from(keys).map_err(ParseError::InvalidFormat))
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
/// An error returned when a vector of keys fails to convert to a format.
pub enum TryFromKeyVectorError {
    /// The input is empty.
    Empty,
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
            Self::Empty => f.write_str("empty input"),
            Self::InvalidGenotypeKeyPosition => f.write_str("invalid genotype key position"),
            Self::DuplicateKey(key) => write!(f, "duplicate key: {}", key),
        }
    }
}

impl TryFrom<Vec<Key>> for Format {
    type Error = TryFromKeyVectorError;

    fn try_from(keys: Vec<Key>) -> Result<Self, Self::Error> {
        if keys.is_empty() {
            return Err(TryFromKeyVectorError::Empty);
        } else if let Some(i) = keys.iter().position(|k| k == &Key::Genotype) {
            if i != 0 {
                return Err(TryFromKeyVectorError::InvalidGenotypeKeyPosition);
            }
        }

        let mut set = HashSet::new();

        for key in &keys {
            if !set.insert(key.clone()) {
                return Err(TryFromKeyVectorError::DuplicateKey(key.clone()));
            }
        }

        Ok(Self(keys))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let format = Format(vec![Key::Genotype]);
        assert_eq!(format.to_string(), "GT");

        let format = Format(vec![
            Key::Genotype,
            Key::ConditionalGenotypeQuality,
            Key::ReadDepth,
            Key::HaplotypeQuality,
        ]);
        assert_eq!(format.to_string(), "GT:GQ:DP:HQ");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("GT".parse(), Ok(Format(vec![Key::Genotype])));
        assert_eq!(
            "GT:GQ".parse(),
            Ok(Format(vec![Key::Genotype, Key::ConditionalGenotypeQuality]))
        );

        assert_eq!("".parse::<Format>(), Err(ParseError::Empty));
        assert!(matches!(
            "GQ:GT".parse::<Format>(),
            Err(ParseError::InvalidFormat(_))
        ));
    }

    #[test]
    fn test_try_from_vec_key_for_format() {
        assert_eq!(
            Format::try_from(vec![Key::Genotype]),
            Ok(Format(vec![Key::Genotype]))
        );

        assert_eq!(
            Format::try_from(vec![Key::Genotype, Key::ConditionalGenotypeQuality]),
            Ok(Format(vec![Key::Genotype, Key::ConditionalGenotypeQuality]))
        );

        assert_eq!(
            Format::try_from(vec![Key::ConditionalGenotypeQuality]),
            Ok(Format(vec![Key::ConditionalGenotypeQuality]))
        );

        assert_eq!(
            Format::try_from(Vec::new()),
            Err(TryFromKeyVectorError::Empty)
        );

        assert_eq!(
            Format::try_from(vec![Key::ConditionalGenotypeQuality, Key::Genotype]),
            Err(TryFromKeyVectorError::InvalidGenotypeKeyPosition)
        );

        assert_eq!(
            Format::try_from(vec![Key::Genotype, Key::Genotype]),
            Err(TryFromKeyVectorError::DuplicateKey(Key::Genotype))
        );
    }
}
