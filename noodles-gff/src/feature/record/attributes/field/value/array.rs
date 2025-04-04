use std::{borrow::Cow, io};

use bstr::BStr;

/// A feature record attributes field array value.
pub trait Array<'r> {
    /// Returns an iterator over values.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Cow<'r, BStr>>> + 'r>;
}
