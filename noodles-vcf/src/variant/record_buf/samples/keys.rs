//! Variant record genotypes keys.

use indexmap::IndexSet;

type Inner = IndexSet<String>;

/// A variant record samples keys buffer.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Keys(Inner);

impl AsRef<Inner> for Keys {
    fn as_ref(&self) -> &Inner {
        &self.0
    }
}

impl AsMut<Inner> for Keys {
    fn as_mut(&mut self) -> &mut Inner {
        &mut self.0
    }
}

impl Extend<String> for Keys {
    fn extend<T: IntoIterator<Item = String>>(&mut self, iter: T) {
        self.0.extend(iter);
    }
}

impl FromIterator<String> for Keys {
    fn from_iter<T: IntoIterator<Item = String>>(iter: T) -> Self {
        let mut keys = Self::default();
        keys.extend(iter);
        keys
    }
}
