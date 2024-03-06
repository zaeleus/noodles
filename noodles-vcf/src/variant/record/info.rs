//! Variant record info fields.

pub mod field;

use std::io;

use self::field::Value;
use crate::Header;

/// Variant record info fields.
pub trait Info {
    /// Returns whether there are any fields.
    fn is_empty(&self) -> bool;

    /// Returns the number of fields.
    fn len(&self) -> usize;

    /// Returns the value of the given key.
    fn get<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
        key: &str,
    ) -> Option<io::Result<Option<Value<'a>>>>;

    /// Returns an iterator over fields.
    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<(&'a str, Option<Value<'a>>)>> + 'a>;
}

impl Info for Box<dyn Info + '_> {
    fn is_empty(&self) -> bool {
        (**self).is_empty()
    }

    fn len(&self) -> usize {
        (**self).len()
    }

    fn get<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
        key: &str,
    ) -> Option<io::Result<Option<Value<'a>>>> {
        (**self).get(header, key)
    }

    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<(&'a str, Option<Value<'a>>)>> + 'a> {
        (**self).iter(header)
    }
}
