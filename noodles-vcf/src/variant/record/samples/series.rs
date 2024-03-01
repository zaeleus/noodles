//! Variant record samples series.

pub mod value;

use std::io;

pub use self::value::Value;
use crate::Header;

/// A variant record samples series.
pub trait Series {
    /// Returns the name.
    fn name<'a, 'h: 'a>(&'a self, header: &'h Header) -> io::Result<&'a str>;

    /// Return the value at the given index.
    fn get<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
        i: usize,
    ) -> Option<Option<io::Result<Value<'a>>>>;

    /// Returns an iterator over values.
    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<Option<Value<'a>>>> + 'a>;
}
