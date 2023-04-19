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

    pub(crate) fn add(&mut self, value: Value) {
        match (self, value) {
            (Self::Unstructured(list), Value::String(s)) => list.push(s),
            (Self::Structured(map), Value::Map(id, m)) => {
                map.insert(id, m);
            }
            (_, _) => panic!("invalid"),
        }
    }
}
