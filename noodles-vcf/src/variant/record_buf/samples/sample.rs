//! VCF record genotype sample.

pub mod value;

use std::{hash::Hash, io};

pub use self::value::Value;
use super::Keys;
use crate::Header;

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
            .as_ref()
            .get_index_of(key)
            .and_then(|i| self.values.get(i).map(|value| value.as_ref()))
    }
}

impl crate::variant::record::samples::Sample for Sample<'_> {
    fn get<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
        key: &str,
    ) -> Option<io::Result<Option<crate::variant::record::samples::series::Value<'a>>>> {
        self.keys
            .as_ref()
            .get_index_of(key)
            .and_then(|i| self.get_index(header, i))
    }

    fn get_index<'a, 'h: 'a>(
        &'a self,
        _: &'h Header,
        i: usize,
    ) -> Option<io::Result<Option<crate::variant::record::samples::series::Value<'a>>>> {
        self.values
            .get(i)
            .map(|value| Ok(value.as_ref().map(|v| v.into())))
    }

    fn iter<'a, 'h: 'a>(
        &'a self,
        _: &'h Header,
    ) -> Box<
        dyn Iterator<
                Item = io::Result<(
                    &'a str,
                    Option<crate::variant::record::samples::series::Value<'a>>,
                )>,
            > + 'a,
    > {
        Box::new(
            self.keys
                .as_ref()
                .iter()
                .zip(self.values)
                .map(|(key, value)| Ok((key.as_ref(), value.as_ref().map(|v| v.into())))),
        )
    }
}
