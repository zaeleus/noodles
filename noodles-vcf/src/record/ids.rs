use std::{fmt, iter};

use memchr::memchr_iter;

use crate::variant::record::Ids as _;

const DELIMITER: char = ';';

/// VCF record IDs.
#[derive(Eq, PartialEq)]
pub struct Ids<'a>(&'a str);

impl<'a> Ids<'a> {
    pub(super) fn new(buf: &'a str) -> Self {
        Self(buf)
    }
}

impl AsRef<str> for Ids<'_> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

impl fmt::Debug for Ids<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl crate::variant::record::Ids for Ids<'_> {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        count(self.0)
    }

    fn iter(&self) -> Box<dyn Iterator<Item = &str> + '_> {
        if self.is_empty() {
            return Box::new(iter::empty());
        }

        Box::new(self.0.split(DELIMITER))
    }
}

fn count(s: &str) -> usize {
    if s.is_empty() {
        0
    } else {
        let n = memchr_iter(DELIMITER as u8, s.as_bytes()).count();
        n + 1
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
    fn test_len() {
        assert_eq!(Ids::new("").len(), 0);
        assert_eq!(Ids::new("nd0").len(), 1);
        assert_eq!(Ids::new("nd0;nd1").len(), 2);
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
