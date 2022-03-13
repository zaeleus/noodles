mod builder;
mod key;

pub use self::{builder::Builder, key::Key};

use std::ops::Deref;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TagIdsDictionary(Vec<Vec<Key>>);

impl Deref for TagIdsDictionary {
    type Target = [Vec<Key>];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<Vec<Vec<Key>>> for TagIdsDictionary {
    fn from(dictionary: Vec<Vec<Key>>) -> Self {
        Self(dictionary)
    }
}
