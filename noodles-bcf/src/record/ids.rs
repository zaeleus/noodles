use std::iter;

/// BCF record IDs.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Ids<'a>(&'a [u8]);

impl<'a> Ids<'a> {
    pub(super) fn new(src: &'a [u8]) -> Self {
        Self(src)
    }

    /// Returns whether there are any IDs.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of IDs.
    pub fn len(&self) -> usize {
        self.iter().count()
    }

    /// Returns an iterator over IDs.
    pub fn iter(&self) -> Box<dyn Iterator<Item = &[u8]> + '_> {
        const DELIMITER: u8 = b';';

        if self.is_empty() {
            Box::new(iter::empty())
        } else {
            Box::new(self.0.split(|&b| b == DELIMITER))
        }
    }
}

impl<'a> AsRef<[u8]> for Ids<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}
