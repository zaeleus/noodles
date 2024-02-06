//! VCF record IDs.

use std::{fmt, ops::Deref, ops::DerefMut};

use indexmap::IndexSet;

const DELIMITER: char = ';';

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

impl fmt::Display for Ids {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, id) in self.iter().enumerate() {
            if i > 0 {
                write!(f, "{DELIMITER}")?;
            }

            f.write_str(id)?;
        }

        Ok(())
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert!(Ids::default().to_string().is_empty());

        let id0 = String::from("nd0");
        let id1 = String::from("nd1");

        let ids: Ids = [id0.clone()].into_iter().collect();
        assert_eq!(ids.to_string(), "nd0");

        let ids: Ids = [id0, id1].into_iter().collect();
        assert_eq!(ids.to_string(), "nd0;nd1");
    }
}
