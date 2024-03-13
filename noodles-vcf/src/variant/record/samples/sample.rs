use std::io;

use super::series::Value;
use crate::Header;

/// Variant record samples sample.
pub trait Sample {
    /// Returns the value with the given key.
    fn get<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
        key: &str,
    ) -> Option<io::Result<Option<Value<'a>>>>;

    /// Returns the value at the given index.
    fn get_index<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
        i: usize,
    ) -> Option<io::Result<Option<Value<'a>>>>;

    /// Returns an iterator over fields.
    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<(&'a str, Option<Value<'a>>)>> + 'a>;
}

impl Sample for Box<dyn Sample + '_> {
    fn get<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
        key: &str,
    ) -> Option<io::Result<Option<Value<'a>>>> {
        (**self).get(header, key)
    }

    fn get_index<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
        i: usize,
    ) -> Option<io::Result<Option<Value<'a>>>> {
        (**self).get_index(header, i)
    }

    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<(&'a str, Option<Value<'a>>)>> + 'a> {
        (**self).iter(header)
    }
}
