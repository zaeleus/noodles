use std::iter;

/// Raw VCF record IDs.
#[derive(Debug, Eq, PartialEq)]
pub struct Ids<'a>(&'a str);

impl<'a> Ids<'a> {
    pub(super) fn new(buf: &'a str) -> Self {
        Self(buf)
    }

    /// Returns whether there are any IDs.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns an iterator over all IDs.
    pub fn iter(&self) -> Box<dyn Iterator<Item = &str> + '_> {
        const DELIMITER: char = ';';

        if self.is_empty() {
            return Box::new(iter::empty());
        }

        Box::new(self.0.split(DELIMITER))
    }
}

impl<'a> AsRef<str> for Ids<'a> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        assert!(Ids::new("").is_empty());
        assert!(!Ids::new("nd0").is_empty());
        assert!(!Ids::new("nd0;nd1").is_empty());
    }

    #[test]
    fn test_iter() {
        let ids = Ids::new("");
        assert!(ids.iter().next().is_none());

        let ids = Ids::new("nd0;nd1");
        let actual: Vec<_> = ids.iter().collect();
        let expected = ["nd0", "nd1"];
        assert_eq!(actual, expected);
    }
}
