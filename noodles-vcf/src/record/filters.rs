use std::{fmt, io};

use crate::Header;

/// VCF record filters.
#[derive(Eq, PartialEq)]
pub struct Filters<'r>(&'r str);

impl<'r> Filters<'r> {
    pub(super) fn new(buf: &'r str) -> Self {
        Self(buf)
    }
}

impl AsRef<str> for Filters<'_> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

impl crate::variant::record::Filters for Filters<'_> {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        let header = Header::default();
        self.iter(&header).count()
    }

    fn iter<'a, 'h: 'a>(
        &'a self,
        _: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<&'a str>> + 'a> {
        iter(self.0)
    }
}

impl fmt::Debug for Filters<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(iter(self.0)).finish()
    }
}

fn iter(s: &str) -> Box<dyn Iterator<Item = io::Result<&str>> + '_> {
    const DELIMITER: char = ';';

    if s.is_empty() {
        Box::new(std::iter::empty())
    } else {
        Box::new(s.split(DELIMITER).map(Ok))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record::Filters as _;

    #[test]
    fn test_is_empty() {
        assert!(Filters::new("").is_empty());
        assert!(!Filters::new("PASS").is_empty());
        assert!(!Filters::new("q10;s50").is_empty());
    }

    #[test]
    fn test_len() {
        assert_eq!(Filters::new("").len(), 0);
        assert_eq!(Filters::new("PASS").len(), 1);
        assert_eq!(Filters::new("q10;s50").len(), 2);
    }

    #[test]
    fn test_iter() -> io::Result<()> {
        let header = Header::default();

        let filters = Filters::new("");
        assert!(filters.iter(&header).next().is_none());

        let filters = Filters::new("PASS");
        let actual: Vec<_> = filters.iter(&header).collect::<io::Result<_>>()?;
        assert_eq!(actual, ["PASS"]);

        let filters = Filters::new("q10;s50");
        let actual: Vec<_> = filters.iter(&header).collect::<io::Result<_>>()?;
        assert_eq!(actual, ["q10", "s50"]);

        Ok(())
    }
}
