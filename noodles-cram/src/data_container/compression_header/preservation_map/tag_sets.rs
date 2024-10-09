//! CRAM data container preservation map tag IDs dictionary

mod builder;
mod key;

pub(crate) use self::builder::Builder;
pub use self::key::Key;

use std::ops::Deref;

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct TagSets(Vec<Vec<Key>>);

impl Deref for TagSets {
    type Target = [Vec<Key>];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<Vec<Vec<Key>>> for TagSets {
    fn from(dictionary: Vec<Vec<Key>>) -> Self {
        Self(dictionary)
    }
}
