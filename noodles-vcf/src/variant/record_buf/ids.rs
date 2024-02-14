//! VCF record IDs.

use indexmap::IndexSet;

/// VCF record IDs (`ID`).
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Ids(IndexSet<String>);

impl AsRef<IndexSet<String>> for Ids {
    fn as_ref(&self) -> &IndexSet<String> {
        &self.0
    }
}

impl AsMut<IndexSet<String>> for Ids {
    fn as_mut(&mut self) -> &mut IndexSet<String> {
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

impl crate::variant::record::Ids for Ids {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = &str> + '_> {
        Box::new(self.0.iter().map(|id| id.as_ref()))
    }
}
