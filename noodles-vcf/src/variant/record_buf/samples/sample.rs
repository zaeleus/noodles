//! VCF record genotype sample.

pub mod value;

pub use self::value::Value;

use std::{error, fmt, hash::Hash};

use super::Keys;

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
