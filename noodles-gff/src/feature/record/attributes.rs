//! Feature record attributes.

pub mod field;

use std::{borrow::Cow, io};

use bstr::BStr;

use self::field::Value;

/// Feature record attributes.
pub trait Attributes {
    /// Returns whether there are any attributes.
    fn is_empty(&self) -> bool;

    /// Returns the value of the given tag.
    fn get(&self, tag: &[u8]) -> Option<io::Result<Value<'_>>>;

    /// Returns an iterator over all tag-value pairs.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Cow<'_, BStr>, Value<'_>)>> + '_>;
}
