mod field;

use std::{io, iter};

use self::field::parse_field;
use crate::{variant::record::info::field::Value, Header};

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
    pub fn iter<'h: 'a>(
        &'a self,
        header: &'h Header,
    ) -> impl Iterator<Item = io::Result<(&'a str, Option<Value<'a>>)>> + 'a {
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
        assert!(Info::new("").is_empty());
        assert!(!Info::new("NS=2;DP=.").is_empty());
    }

    #[test]
    fn test_iter() {
        use crate::record::info::field::key;

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
