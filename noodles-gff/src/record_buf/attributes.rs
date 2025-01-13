//! GFF record attributes.

pub mod field;

use std::ops::Deref;

use indexmap::IndexMap;

use self::field::{Tag, Value};

/// GFF record attributes.
///
/// Attributes are extra data attached to a GFF record. They are represented as a typed map, where
/// each key ([`Tag`]) is associated with a typed [`Value`].
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Attributes(IndexMap<Tag, Value>);

impl Deref for Attributes {
    type Target = IndexMap<Tag, Value>;

    fn deref(&self) -> &Self::Target {
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
