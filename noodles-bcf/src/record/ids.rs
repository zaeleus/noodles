use std::iter;

/// BCF record IDs.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Ids<'a>(&'a [u8]);

impl<'a> Ids<'a> {
    pub(super) fn new(src: &'a [u8]) -> Self {
        Self(src)
    }

    /// Returns whether there are any IDs.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of IDs.
    pub fn len(&self) -> usize {
        self.iter().count()
    }

    /// Returns an iterator over IDs.
    pub fn iter(&self) -> Box<dyn Iterator<Item = &[u8]> + '_> {
        const DELIMITER: u8 = b';';

        if self.is_empty() {
            Box::new(iter::empty())
        } else {
            Box::new(self.0.split(|&b| b == DELIMITER))
        }
    }
}

impl<'a> AsRef<[u8]> for Ids<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        assert!(Ids::new(b"").is_empty());
        assert!(!Ids::new(b"nd0").is_empty());
        assert!(!Ids::new(b"nd0;nd1").is_empty());
    }

    #[test]
    fn test_len() {
        assert_eq!(Ids::new(b"").len(), 0);
        assert_eq!(Ids::new(b"nd0").len(), 1);
        assert_eq!(Ids::new(b"nd0;nd1").len(), 2);
    }

    #[test]
    fn test_iter() {
        let ids = Ids::new(b"");
        assert!(ids.iter().next().is_none());

        let ids = Ids::new(b"nd0;nd1");
        let actual: Vec<_> = ids.iter().collect();
        let expected = [b"nd0", b"nd1"];
        assert_eq!(actual, expected);
    }
}
