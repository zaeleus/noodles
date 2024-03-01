use std::io;

use super::series::Value;
use crate::Header;

/// Variant record samples sample.
pub trait Sample {
    /// Returns an iterator over fields.
    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<(&'a str, Option<Value<'a>>)>> + 'a>;
}

impl Sample for Box<dyn Sample + '_> {
    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<(&'a str, Option<Value<'a>>)>> + 'a> {
        (**self).iter(header)
    }
}
