use std::io;

use super::Keys;
use crate::variant::record::samples::series::Value;

/// A VCF record samples sample.
#[derive(Debug, Eq, PartialEq)]
pub struct Sample<'a> {
    src: &'a str,
    keys: Keys<'a>,
}

impl<'a> Sample<'a> {
    pub(super) fn new(src: &'a str, keys: Keys<'a>) -> Self {
        Self { src, keys }
    }

    pub fn values(&self) -> impl Iterator<Item = Option<io::Result<Value<'_>>>> + '_ {
        const MISSING: &str = ".";
        const DELIMITER: char = ':';

        self.src.split(DELIMITER).map(|s| match s {
            MISSING => None,
            _ => Some(parse_value(s)),
        })
    }

    /// Returns an iterator over fields.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<(&str, Option<Value<'_>>)>> + '_ {
        self.keys
            .iter()
            .zip(self.values())
            .map(|(key, value)| value.transpose().map(|v| (key, v)))
    }
}

impl<'a> AsRef<str> for Sample<'a> {
    fn as_ref(&self) -> &str {
        self.src
    }
}

fn parse_value(src: &str) -> io::Result<Value<'_>> {
    Ok(Value::String(src))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_values() {
        let keys = Keys::new("GT:GQ");
        let sample = Sample::new("0|0:.", keys);
        let mut iter = sample.values();

        assert!(matches!(
            iter.next(),
            Some(Some(Ok(Value::String(s)))) if s == "0|0"
        ));

        assert!(matches!(iter.next(), Some(None)));

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_iter() {
        let keys = Keys::new("GT:GQ");
        let sample = Sample::new("0|0:.", keys);
        let mut iter = sample.iter();

        assert!(matches!(
            iter.next(),
            Some(Ok(("GT", Some(Value::String(s))))) if s == "0|0"
        ));

        assert!(matches!(iter.next(), Some(Ok(("GQ", None)))));

        assert!(iter.next().is_none());
    }
}
