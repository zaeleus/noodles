use std::iter;

use crate::record::MISSING_FIELD;

/// Raw VCF record info.
#[derive(Debug, Eq, PartialEq)]
pub struct Info<'a>(&'a str);

impl<'a> Info<'a> {
    pub(super) fn new(buf: &'a str) -> Self {
        Self(buf)
    }

    /// Returns whether there are any info fields.
    pub fn is_empty(&self) -> bool {
        self.0 == MISSING_FIELD
    }

    /// Returns an iterator over all fields.
    pub fn iter(&self) -> Box<dyn Iterator<Item = (&str, Option<&str>)> + '_> {
        const DELIMITER: char = ';';
        const FIELD_SEPARATOR: char = '=';

        if self.is_empty() {
            return Box::new(iter::empty());
        }

        Box::new(self.0.split(DELIMITER).map(|raw_field| {
            let mut components = raw_field.split(FIELD_SEPARATOR);
            let key = components.next().unwrap();
            let value = components.next();
            (key, value)
        }))
    }
}

impl<'a> AsRef<str> for Info<'a> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        assert!(Info::new(".").is_empty());
        assert!(!Info::new("NS=2;DP=.").is_empty());
    }

    #[test]
    fn test_iter() {
        let info = Info::new(".");
        assert!(info.iter().next().is_none());

        let info = Info::new("NS=2;DP=.");
        let actual: Vec<_> = info.iter().collect();
        let expected = [("NS", Some("2")), ("DP", None)];
        assert_eq!(actual, expected);
    }
}
