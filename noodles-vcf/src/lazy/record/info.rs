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
        self.0.is_empty()
    }

    /// Returns an iterator over all fields.
    pub fn iter(&self) -> impl Iterator<Item = (&str, Option<&str>)> + '_ {
        let mut src = self.0;

        iter::from_fn(move || {
            if src.is_empty() {
                None
            } else {
                Some(parse_field(&mut src))
            }
        })
    }
}

impl<'a> AsRef<str> for Info<'a> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

fn parse_field<'a>(src: &mut &'a str) -> (&'a str, Option<&'a str>) {
    const DELIMITER: u8 = b';';
    const FIELD_SEPARATOR: char = '=';

    let raw_field = match src.as_bytes().iter().position(|&b| b == DELIMITER) {
        Some(i) => {
            let (buf, rest) = src.split_at(i);
            *src = &rest[1..];
            buf
        }
        None => {
            let (buf, rest) = src.split_at(src.len());
            *src = rest;
            buf
        }
    };

    let mut components = raw_field.split(FIELD_SEPARATOR);
    let key = components.next().unwrap();
    let value = components.next().and_then(|s| match s {
        MISSING_FIELD => None,
        _ => Some(s),
    });

    (key, value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        assert!(Info::new("").is_empty());
        assert!(!Info::new("NS=2;DP=.").is_empty());
    }

    #[test]
    fn test_iter() {
        let info = Info::new("");
        assert!(info.iter().next().is_none());

        let info = Info::new("NS=2;DP=.");
        let actual: Vec<_> = info.iter().collect();
        let expected = [("NS", Some("2")), ("DP", None)];
        assert_eq!(actual, expected);
    }
}
