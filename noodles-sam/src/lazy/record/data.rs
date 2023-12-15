use std::{io, iter};

mod field;

use self::field::{parse_field, Tag};
use crate::alignment::record::data::field::Value;

/// Raw SAM record data.
pub struct Data<'a>(&'a [u8]);

impl<'a> Data<'a> {
    pub(super) fn new(buf: &'a [u8]) -> Self {
        Self(buf)
    }

    /// Returns whether there are any data fields.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
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

impl<'a> AsRef<[u8]> for Data<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iter() -> io::Result<()> {
        let data = Data::new(b"");
        assert!(data.iter().next().is_none());

        let data = Data::new(b"NH:i:1");
        let actual: Vec<_> = data.iter().collect::<io::Result<_>>()?;

        assert_eq!(actual.len(), 1);

        let (actual_tag, actual_value) = &actual[0];
        assert_eq!(actual_tag, &[b'N', b'H']);
        assert!(matches!(actual_value, Value::Int32(1)));

        Ok(())
    }
}
