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
        self.values(header).nth(i)
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

        if self.as_ref().is_empty() {
            Box::new(iter::empty())
        } else {
            Box::new(
                self.keys
                    .iter()
                    .zip(self.src.split(DELIMITER))
                    .map(|(key, s)| parse_value(s, header, key).map(|value| (key, value))),
            )
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
        self.iter(header)
            .nth(i)
            .map(|result| result.map(|(_, value)| value))
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
}
