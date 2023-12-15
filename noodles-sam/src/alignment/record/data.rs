//! Alignment record data.

#![allow(dead_code)]

pub mod field;

use std::io;

use self::field::Value;

/// Alignment record data.
pub trait Data {
    /// Returns whether there are any fields.
    fn is_empty(&self) -> bool;

    /// Returns the value for the given tag.
    fn get(&self, tag: &[u8; 2]) -> Option<io::Result<Value<'_>>>;

    /// Returns an iterator over fields.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<([u8; 2], Value<'_>)>> + '_>;
}
