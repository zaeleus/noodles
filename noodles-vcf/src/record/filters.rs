use std::iter;

/// Raw VCF record filters.
#[derive(Debug, Eq, PartialEq)]
pub struct Filters<'a>(&'a str);

impl<'a> Filters<'a> {
    pub(super) fn new(buf: &'a str) -> Self {
        Self(buf)
    }
}

impl<'a> AsRef<str> for Filters<'a> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

impl<'a> crate::variant::record::Filters for Filters<'a> {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        self.iter().count()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = &str> + '_> {
        const DELIMITER: char = ';';

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
    fn test_iter() {
        let filters = Filters::new("");
        assert!(filters.iter().next().is_none());

        let filters = Filters::new("PASS");
        let actual: Vec<_> = filters.iter().collect();
        assert_eq!(actual, ["PASS"]);

        let filters = Filters::new("q10;s50");
        let actual: Vec<_> = filters.iter().collect();
        assert_eq!(actual, ["q10", "s50"]);
    }
}
