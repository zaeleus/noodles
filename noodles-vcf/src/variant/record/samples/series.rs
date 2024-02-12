//! Variant record samples series.

pub mod value;

use std::io;

pub use self::value::Value;

/// A variant record samples series.
pub trait Series {
    /// Returns the name.
    fn name(&self) -> &str;

    /// Returns an iterator over values.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Value<'_>>> + '_>;
}
