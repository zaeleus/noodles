use std::{fmt, io, iter};

use crate::variant::record::AlternateBases as _;

/// VCF record alternate bases.
pub struct AlternateBases<'a>(&'a str);

const DELIMITER: char = ',';

impl<'a> AlternateBases<'a> {
    pub(super) fn new(src: &'a str) -> Self {
        Self(src)
    }
}

impl fmt::Debug for AlternateBases<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl AsRef<str> for AlternateBases<'_> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

impl crate::variant::record::AlternateBases for AlternateBases<'_> {
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

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<&str>> + '_> {
        if self.is_empty() {
            Box::new(iter::empty())
        } else {
            Box::new(self.0.split(DELIMITER).map(Ok))
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
    fn test_iter() -> io::Result<()> {
        let alternate_bases = AlternateBases::new("");
        assert!(alternate_bases.iter().next().is_none());

        let alternate_bases = AlternateBases::new("A");
        let actual: Vec<_> = alternate_bases.iter().collect::<io::Result<_>>()?;
        let expected = ["A"];
        assert_eq!(actual, expected);

        let alternate_bases = AlternateBases::new("A,C");
        let actual: Vec<_> = alternate_bases.iter().collect::<io::Result<_>>()?;
        let expected = ["A", "C"];
        assert_eq!(actual, expected);

        Ok(())
    }
}
