//! VCF record IDs.

use std::{ops::Deref, ops::DerefMut};

use indexmap::IndexSet;

/// VCF record IDs (`ID`).
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Ids(IndexSet<String>);

impl Deref for Ids {
    type Target = IndexSet<String>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Ids {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Extend<String> for Ids {
    fn extend<T: IntoIterator<Item = String>>(&mut self, iter: T) {
        self.0.extend(iter);
    }
}

impl FromIterator<String> for Ids {
    fn from_iter<T: IntoIterator<Item = String>>(iter: T) -> Self {
        let mut ids = Self::default();
        ids.extend(iter);
        ids
    }
}
