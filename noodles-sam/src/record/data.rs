//! SAM record data.

pub mod field;

use std::{borrow::Borrow, fmt, io, iter};

use self::field::parse_field;
use crate::alignment::record::data::field::{Tag, Value};

/// SAM record data.
pub struct Data<'a>(&'a [u8]);

impl<'a> Data<'a> {
    pub(super) fn new(buf: &'a [u8]) -> Self {
        Self(buf)
    }

    /// Returns whether there are any data fields.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the value for the given tag.
    pub fn get<K>(&self, tag: &K) -> Option<io::Result<Value<'_>>>
    where
        K: Borrow<[u8; 2]>,
    {
        for result in self.iter() {
            match result {
                Ok((t, value)) => {
                    if &t == tag.borrow() {
                        return Some(Ok(value));
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        None
    }

    /// Returns an iterator over all tag-value pairs.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<(Tag, Value<'_>)>> + '_ {
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

impl fmt::Debug for Data<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut formatter = f.debug_map();

        for result in self.iter() {
            let (tag, value) = result.map_err(|_| fmt::Error)?;
            formatter.entry(&tag, &value);
        }

        formatter.finish()
    }
}

impl crate::alignment::record::Data for Data<'_> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn get(&self, tag: &Tag) -> Option<io::Result<Value<'_>>> {
        self.get(tag)
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Tag, Value<'_>)>> + '_> {
        Box::new(self.iter())
    }
}

impl AsRef<[u8]> for Data<'_> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        let data = Data::new(b"");
        assert!(data.is_empty());

        let data = Data::new(b"NH:i:1");
        assert!(!data.is_empty());
    }

    #[test]
    fn test_get() -> io::Result<()> {
        let data = Data::new(b"");
        assert!(data.get(&Tag::ALIGNMENT_HIT_COUNT).is_none());
        assert!(data.get(&Tag::COMMENT).is_none());

        let data = Data::new(b"NH:i:1");
        assert!(matches!(
            data.get(&Tag::ALIGNMENT_HIT_COUNT).transpose()?,
            Some(Value::Int32(1))
        ));
        assert!(data.get(&Tag::COMMENT).is_none());

        let data = Data::new(b"NH:i:1\tCO:Z:ndls");
        assert!(matches!(
            data.get(&Tag::ALIGNMENT_HIT_COUNT).transpose()?,
            Some(Value::Int32(1))
        ));
        assert!(matches!(
            data.get(&Tag::COMMENT).transpose()?,
            Some(Value::String(s)) if s == b"ndls"
        ));

        Ok(())
    }

    #[test]
    fn test_iter() -> io::Result<()> {
        let data = Data::new(b"");
        assert!(data.iter().next().is_none());

        let data = Data::new(b"NH:i:1");
        let actual: Vec<_> = data.iter().collect::<io::Result<_>>()?;

        assert_eq!(actual.len(), 1);

        let (actual_tag, actual_value) = &actual[0];
        assert_eq!(actual_tag, &Tag::ALIGNMENT_HIT_COUNT);
        assert!(matches!(actual_value, Value::Int32(1)));

        Ok(())
    }
}
