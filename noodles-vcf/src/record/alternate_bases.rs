use std::{fmt, iter};

use crate::variant::record::AlternateBases as _;

/// VCF record alternate bases.
pub struct AlternateBases<'a>(&'a str);

const DELIMITER: char = ',';

impl<'a> AlternateBases<'a> {
    pub(super) fn new(src: &'a str) -> Self {
        Self(src)
    }
}

impl<'a> fmt::Debug for AlternateBases<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl<'a> AsRef<str> for AlternateBases<'a> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

impl<'a> crate::variant::record::AlternateBases for AlternateBases<'a> {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        if self.is_empty() {
            0
        } else {
            self.iter().count()
        }
    }

    fn iter(&self) -> Box<dyn Iterator<Item = &str> + '_> {
        if self.is_empty() {
            Box::new(iter::empty())
        } else {
            Box::new(self.0.split(DELIMITER))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        assert!(AlternateBases::new("").is_empty());
        assert!(!AlternateBases::new("A").is_empty());
        assert!(!AlternateBases::new("A,C").is_empty());
    }

    #[test]
    fn test_len() {
        assert_eq!(AlternateBases::new("").len(), 0);
        assert_eq!(AlternateBases::new("A").len(), 1);
        assert_eq!(AlternateBases::new("A,C").len(), 2);
    }

    #[test]
    fn test_iter() {
        let alternate_bases = AlternateBases::new("");
        assert!(alternate_bases.iter().next().is_none());

        let alternate_bases = AlternateBases::new("A");
        let mut iter = alternate_bases.iter();
        assert_eq!(iter.next(), Some("A"));
        assert!(iter.next().is_none());

        let alternate_bases = AlternateBases::new("A,C");
        let mut iter = alternate_bases.iter();
        assert_eq!(iter.next(), Some("A"));
        assert_eq!(iter.next(), Some("C"));
        assert!(iter.next().is_none());
    }
}
