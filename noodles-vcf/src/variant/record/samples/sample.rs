use std::io;

use super::series::Value;

/// Variant record samples sample.
pub trait Sample {
    /// Returns an iterator over fields.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(&str, Option<Value<'_>>)>> + '_>;
}
