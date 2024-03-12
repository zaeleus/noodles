use std::{io, str};

/// Variant record reference bases.
pub trait ReferenceBases {
    /// Returns whether there are any reference bases.
    fn is_empty(&self) -> bool;

    /// Returns the number of reference bases.
    fn len(&self) -> usize;

    /// Returns an iterator over reference bases.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u8>> + '_>;
}

impl ReferenceBases for &str {
    fn is_empty(&self) -> bool {
        str::is_empty(self)
    }

    fn len(&self) -> usize {
        str::len(self)
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u8>> + '_> {
        Box::new(self.as_bytes().iter().copied().map(Ok))
    }
}

impl ReferenceBases for Box<dyn ReferenceBases + '_> {
    fn is_empty(&self) -> bool {
        (**self).is_empty()
    }

    fn len(&self) -> usize {
        (**self).len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u8>> + '_> {
        (**self).iter()
    }
}
