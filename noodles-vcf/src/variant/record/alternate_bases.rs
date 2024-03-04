use std::io;

/// Variant record alternate bases.
pub trait AlternateBases {
    /// Returns whether there are any alternate bases.
    fn is_empty(&self) -> bool;

    /// Returns the number of alternate bases.
    fn len(&self) -> usize;

    /// Returns an iterator over alternate bases.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<&str>> + '_>;
}

impl AlternateBases for Box<dyn AlternateBases + '_> {
    fn is_empty(&self) -> bool {
        (**self).is_empty()
    }

    fn len(&self) -> usize {
        (**self).len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<&str>> + '_> {
        (**self).iter()
    }
}
