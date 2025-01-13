//! GFF record attributes.

pub mod field;

use indexmap::IndexMap;

use self::field::{Tag, Value};

/// GFF record attributes.
///
/// Attributes are extra data attached to a GFF record. They are represented as a typed map, where
/// each key ([`Tag`]) is associated with a typed [`Value`].
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Attributes(IndexMap<Tag, Value>);

impl Attributes {
    /// Returns whether there are any entries.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of entries.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns the value at the given tag.
    pub fn get(&self, tag: &str) -> Option<&Value> {
        self.0.get(tag)
    }
}

impl AsRef<IndexMap<Tag, Value>> for Attributes {
    fn as_ref(&self) -> &IndexMap<Tag, Value> {
        &self.0
    }
}

impl AsMut<IndexMap<Tag, Value>> for Attributes {
    fn as_mut(&mut self) -> &mut IndexMap<Tag, Value> {
        &mut self.0
    }
}

impl Extend<(Tag, Value)> for Attributes {
    fn extend<T: IntoIterator<Item = (Tag, Value)>>(&mut self, iter: T) {
        self.0.extend(iter);
    }
}

impl FromIterator<(Tag, Value)> for Attributes {
    fn from_iter<T: IntoIterator<Item = (Tag, Value)>>(iter: T) -> Self {
        let mut attributes = Self::default();
        attributes.extend(iter);
        attributes
    }
}
