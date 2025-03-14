//! Feature record attributes.

pub mod field;

use std::{borrow::Cow, io};

use self::field::Value;

/// Feature record attributes.
pub trait Attributes {
    /// Returns whether there are any attributes.
    fn is_empty(&self) -> bool;

    /// Returns the value fo the given tag.
    fn get(&self, tag: &str) -> Option<io::Result<Value<'_>>>;

    /// Returns an iterator over all tag-value pairs.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Cow<'_, str>, Value<'_>)>> + '_>;
}
