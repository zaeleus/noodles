//! VCF record genotype sample.

pub mod value;

use std::{hash::Hash, io};

pub use self::value::Value;
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

impl<'g> crate::variant::record::samples::Sample for Sample<'g> {
    fn iter(
        &self,
    ) -> Box<
        dyn Iterator<
                Item = io::Result<(
                    &str,
                    Option<crate::variant::record::samples::series::Value<'_>>,
                )>,
            > + '_,
    > {
        Box::new(
            self.keys
                .iter()
                .zip(self.values)
                .map(|(key, value)| Ok((key.as_ref(), value.as_ref().map(|v| v.into())))),
        )
    }
}
