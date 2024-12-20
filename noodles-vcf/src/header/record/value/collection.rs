//! VCF header record value collection.

use std::{error, fmt};

use indexmap::IndexMap;

use super::Value;
use crate::header::record::value::{map, Map};

/// A collection of VCF header record other values.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Collection {
    /// A list of unstructured strings.
    Unstructured(Vec<String>),
    /// A map of structured maps.
    Structured(IndexMap<String, Map<map::Other>>),
}

impl Collection {
    /// Returns whether the collection has values.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the number of values in the collection.
    pub fn len(&self) -> usize {
        match self {
            Self::Unstructured(list) => list.len(),
            Self::Structured(map) => map.len(),
        }
    }

    pub(crate) fn add(&mut self, value: Value) -> Result<(), AddError> {
        match (self, value) {
            (Self::Unstructured(list), Value::String(s)) => {
                list.push(s);
                Ok(())
            }
            (Self::Unstructured(_), Value::Map(..)) => Err(AddError::TypeMismatch {
                actual: "structured",
                expected: "unstructured",
            }),
            (Self::Structured(map), Value::Map(id, m)) => try_insert(map, id, m),
            (Self::Structured(_), Value::String(_)) => Err(AddError::TypeMismatch {
                actual: "unstructured",
                expected: "structured",
            }),
        }
    }
}

fn try_insert(
    map: &mut IndexMap<String, Map<map::Other>>,
    id: String,
    m: Map<map::Other>,
) -> Result<(), AddError> {
    use indexmap::map::Entry;

    match map.entry(id) {
        Entry::Vacant(entry) => {
            entry.insert(m);
            Ok(())
        }
        Entry::Occupied(entry) => Err(AddError::DuplicateId(entry.key().into())),
    }
}

/// An error returned when a value fails be added to a collection.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum AddError {
    /// The value does not match the collection type.
    TypeMismatch {
        /// The value type.
        actual: &'static str,
        /// The collection type.
        expected: &'static str,
    },
    /// An ID is duplicated.
    DuplicateId(String),
}

impl error::Error for AddError {}

impl fmt::Display for AddError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TypeMismatch { actual, expected } => {
                write!(f, "type mismatch: expected {expected}, got {actual}")
            }
            Self::DuplicateId(id) => write!(f, "duplicate ID: {id}"),
        }
    }
}
