//! VCF record genotype sample.

pub mod value;

pub use self::value::Value;

use std::{error, fmt, hash::Hash};

use super::{keys::key, Keys};

/// A VCF record genotype sample.
#[derive(Debug, PartialEq)]
pub struct Sample<'g> {
    keys: &'g Keys,
    values: &'g [Option<Value>],
}

impl<'g> Sample<'g> {
    /// Creates a new genotype sample.
    pub fn new(keys: &'g Keys, values: &'g [Option<Value>]) -> Self {
        Self { keys, values }
    }

    /// Returns the keys.
    pub fn keys(&self) -> &'g Keys {
        self.keys
    }

    /// Returns the values.
    pub fn values(&self) -> &'g [Option<Value>] {
        self.values
    }

    /// Returns a reference to the value with the given key.
    pub fn get<K>(&self, key: &K) -> Option<Option<&'g Value>>
    where
        K: Hash + indexmap::Equivalent<String> + ?Sized,
    {
        self.keys
            .get_index_of(key)
            .and_then(|i| self.values.get(i).map(|value| value.as_ref()))
    }

    /// Returns the VCF record genotypes genotype value.
    ///
    /// This is a convenience method to return a parsed version of the genotype (`GT`) field value.
    pub fn genotype(&self) -> Option<Result<value::Genotype, GenotypeError>> {
        self.get(key::GENOTYPE).map(|value| match value {
            Some(Value::String(s)) => s.parse().map_err(GenotypeError::InvalidValue),
            _ => Err(GenotypeError::InvalidValueType(value.cloned())),
        })
    }
}

/// An error returned when a raw VCF genotype fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// A value is invalid.
    InvalidValue(value::ParseError),
    /// The genotype field value is unexpected.
    ///
    /// There are fewer keys than values.
    UnexpectedValue,
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidValue(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::InvalidValue(_) => f.write_str("invalid value"),
            Self::UnexpectedValue => f.write_str("unexpected value"),
        }
    }
}

/// An error returned when a genotype (`GT`) field value fails to parse.
#[derive(Clone, Debug, PartialEq)]
pub enum GenotypeError {
    /// The genotype field value is invalid.
    InvalidValue(value::genotype::ParseError),
    /// The genotype field value type is invalid.
    ///
    /// The `GT` field value must be a `String`.
    InvalidValueType(Option<Value>),
}

impl error::Error for GenotypeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidValue(e) => Some(e),
            Self::InvalidValueType(_) => None,
        }
    }
}

impl fmt::Display for GenotypeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidValue(_) => f.write_str("invalid value"),
            Self::InvalidValueType(value) => write!(f, "invalid String, got {value:?}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genotype() -> Result<(), crate::record::genotypes::keys::TryFromKeyVectorError> {
        let keys = Keys::try_from(vec![String::from(key::GENOTYPE)])?;

        let values = vec![Some(Value::from("ndls"))];
        let sample = Sample::new(&keys, &values);

        assert!(matches!(
            sample.genotype(),
            Some(Err(GenotypeError::InvalidValue(_)))
        ));

        let values = vec![Some(Value::from(0))];
        let sample = Sample::new(&keys, &values);

        assert!(matches!(
            sample.genotype(),
            Some(Err(GenotypeError::InvalidValueType(_)))
        ));

        Ok(())
    }
}
