use std::{io, iter};

use super::{Keys, series::value::parse_value};
use crate::{Header, variant::record::samples::series::Value};

/// A VCF record samples sample.
#[derive(Debug, Eq, PartialEq)]
pub struct Sample<'s> {
    src: &'s str,
    keys: Keys<'s>,
}

impl<'s> Sample<'s> {
    pub(super) fn new(src: &'s str, keys: Keys<'s>) -> Self {
        Self { src, keys }
    }

    /// Returns the value at the given index.
    pub fn get_index<'h: 's>(
        &self,
        header: &'h Header,
        i: usize,
    ) -> Option<Option<io::Result<Value<'s>>>> {
        if self.src.is_empty() {
            Some(Some(Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "empty sample",
            ))))
        } else {
            self.values(header).nth(i)
        }
    }

    /// Returns an iterator over values.
    pub fn values<'h: 's>(
        &self,
        header: &'h Header,
    ) -> impl Iterator<Item = Option<io::Result<Value<'s>>>> + '_ {
        self.iter(header)
            .map(|result| result.map(|(_, value)| value).transpose())
    }

    /// Returns an iterator over fields.
    pub fn iter<'h: 's>(
        &self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<(&str, Option<Value<'s>>)>> + '_> {
        const DELIMITER: char = ':';
        const MISSING: &str = ".";

        if self.src.is_empty() {
            Box::new(iter::once(Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "empty field",
            ))))
        } else if self.src == MISSING {
            Box::new(iter::empty())
        } else {
            let mut keys = self.keys.iter();
            let mut values = self.src.split(DELIMITER);
            let mut has_unexpected_value = false;

            Box::new(iter::from_fn(move || {
                if has_unexpected_value {
                    return None;
                }

                match (keys.next(), values.next()) {
                    (Some(key), Some(s)) => {
                        Some(parse_value(s, header, key).map(|value| (key, value)))
                    }
                    (None, Some(_)) => {
                        has_unexpected_value = true;

                        Some(Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "unexpected value",
                        )))
                    }
                    _ => None,
                }
            }))
        }
    }
}

impl AsRef<str> for Sample<'_> {
    fn as_ref(&self) -> &str {
        self.src
    }
}

impl crate::variant::record::samples::Sample for Sample<'_> {
    fn get<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
        key: &str,
    ) -> Option<io::Result<Option<Value<'a>>>> {
        for result in self.iter(header) {
            match result {
                Ok((k, v)) => {
                    if k == key {
                        return Some(Ok(v));
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        None
    }

    fn get_index<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
        i: usize,
    ) -> Option<io::Result<Option<Value<'a>>>> {
        self.get_index(header, i).map(|v| v.transpose())
    }

    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<(&'a str, Option<Value<'a>>)>> + 'a> {
        Box::new(self.iter(header))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_index() {
        let header = Header::default();
        let keys = Keys::new("GT:GQ");
        let sample = Sample::new("0|0:.", keys);

        assert!(matches!(
            sample.get_index(&header, 0),
            Some(Some(Ok(Value::Genotype(_))))
        ));
        assert!(matches!(sample.get_index(&header, 1), Some(None)));
        assert!(sample.get_index(&header, 2).is_none());
    }

    #[test]
    fn test_get_index_with_empty_src() {
        let header = Header::default();
        let keys = Keys::new("GT:GQ");
        let sample = Sample::new("", keys);

        assert!(matches!(
            sample.get_index(&header, 0),
            Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
        ));

        assert!(matches!(
            sample.get_index(&header, 1),
            Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[test]
    fn test_values() {
        let header = Header::default();

        let keys = Keys::new("GT:GQ");
        let sample = Sample::new("0|0:.", keys);
        let mut iter = sample.values(&header);

        assert!(matches!(iter.next(), Some(Some(Ok(Value::Genotype(_))))));
        assert!(matches!(iter.next(), Some(None)));
        assert!(iter.next().is_none());
    }

    #[test]
    fn test_iter() {
        let header = Header::default();

        let keys = Keys::new("GT:GQ");
        let sample = Sample::new("0|0:.", keys);
        let mut iter = sample.iter(&header);

        assert!(matches!(
            iter.next(),
            Some(Ok(("GT", Some(Value::Genotype(_)))))
        ));

        assert!(matches!(iter.next(), Some(Ok(("GQ", None)))));

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_iter_with_missing_field() {
        let header = Header::default();

        let keys = Keys::new("GT:GQ");
        let sample = Sample::new(".", keys);
        let mut iter = sample.iter(&header);

        assert!(iter.next().is_none());
    }

    #[test]
    fn test_iter_with_empty_src() {
        let header = Header::default();

        let keys = Keys::new("GT:GQ");
        let sample = Sample::new("", keys);
        let mut iter = sample.iter(&header);

        assert!(matches!(
            iter.next(),
            Some(Err(e)) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[test]
    fn test_iter_with_unexpected_value() {
        let header = Header::default();

        let keys = Keys::new("GT");
        let sample = Sample::new("0|0:13", keys);
        let mut iter = sample.iter(&header);

        assert!(matches!(
            iter.next(),
            Some(Ok(("GT", Some(Value::Genotype(_)))))
        ));

        assert!(matches!(
            iter.next(),
            Some(Err(e)) if e.kind() == io::ErrorKind::InvalidData
        ));

        assert!(iter.next().is_none());
    }
}
