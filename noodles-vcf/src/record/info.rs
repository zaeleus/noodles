mod field;

use std::{io, iter};

use self::field::parse_field;
use crate::{variant::record::info::field::Value, Header};

/// Raw VCF record info.
#[derive(Debug, Eq, PartialEq)]
pub struct Info<'r>(&'r str);

impl<'r> Info<'r> {
    pub(super) fn new(buf: &'r str) -> Self {
        Self(buf)
    }

    /// Returns the value with the given key.
    pub fn get<'h: 'r>(
        &'r self,
        header: &'h Header,
        key: &str,
    ) -> Option<io::Result<Option<Value<'r>>>> {
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

    /// Returns an iterator over all fields.
    pub fn iter<'h: 'r>(
        &'r self,
        header: &'h Header,
    ) -> impl Iterator<Item = io::Result<(&'r str, Option<Value<'r>>)>> + 'r {
        let mut src = self.0;

        iter::from_fn(move || {
            if src.is_empty() {
                None
            } else {
                Some(parse_field(&mut src, header))
            }
        })
    }
}

impl AsRef<str> for Info<'_> {
    fn as_ref(&self) -> &str {
        self.0
    }
}

impl crate::variant::record::Info for Info<'_> {
    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn len(&self) -> usize {
        const DELIMITER: u8 = b';';

        if self.is_empty() {
            0
        } else {
            let n = self
                .0
                .as_bytes()
                .iter()
                .filter(|&&b| b == DELIMITER)
                .count();

            n + 1
        }
    }

    fn get<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
        key: &str,
    ) -> Option<io::Result<Option<Value<'a>>>> {
        self.get(header, key)
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
    use crate::variant::record::Info as _;

    #[test]
    fn test_is_empty() {
        assert!(Info::new("").is_empty());
        assert!(!Info::new("NS=2;DP=.").is_empty());
    }

    #[test]
    fn test_iter() {
        use crate::variant::record::info::field::key;

        let header = Header::default();

        let info = Info::new("");
        assert!(info.iter(&header).next().is_none());

        let info = Info::new("NS=2;DP=.");
        let mut iter = info.iter(&header);

        assert!(matches!(
            iter.next(),
            Some(Ok((key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2)))))
        ));

        assert!(matches!(iter.next(), Some(Ok((key::TOTAL_DEPTH, None)))));
    }
}
