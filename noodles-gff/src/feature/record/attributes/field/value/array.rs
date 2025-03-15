use std::{borrow::Cow, io};

/// A feature record attributes field array value.
pub trait Array<'r> {
    /// Returns an iterator over values.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Cow<'r, str>>> + 'r>;
}
