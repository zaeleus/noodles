use std::ops::Deref;

use crate::record;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TagIdsDictionary(Vec<Vec<record::tag::Key>>);

impl Deref for TagIdsDictionary {
    type Target = [Vec<record::tag::Key>];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<Vec<Vec<record::tag::Key>>> for TagIdsDictionary {
    fn from(dictionary: Vec<Vec<record::tag::Key>>) -> Self {
        Self(dictionary)
    }
}
