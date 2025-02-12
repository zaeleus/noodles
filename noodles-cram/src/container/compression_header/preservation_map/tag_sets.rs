//! CRAM container preservation map tag sets.

mod key;

pub use self::key::Key;

use std::ops::Deref;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
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
