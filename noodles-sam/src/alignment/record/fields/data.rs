use std::io;

use crate::alignment::record::data::field::{Tag, Value};

/// Alignment record data.
pub trait Data {
    /// Returns whether there are any fields.
    fn is_empty(&self) -> bool;

    /// Returns the value for the given tag.
    fn get(&self, tag: &Tag) -> Option<io::Result<Value<'_>>>;

    /// Returns an iterator over fields.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Tag, Value<'_>)>> + '_>;
}

impl Data for Box<dyn Data + '_> {
    fn is_empty(&self) -> bool {
        (**self).is_empty()
    }

    fn get(&self, tag: &Tag) -> Option<io::Result<Value<'_>>> {
        (**self).get(tag)
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Tag, Value<'_>)>> + '_> {
        (**self).iter()
    }
}
