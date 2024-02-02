/// BCF record IDs.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Ids<'a>(&'a [u8]);

const DELIMITER: u8 = b';';

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
        if self.is_empty() {
            0
        } else {
            self.0.iter().filter(|&&b| b == DELIMITER).count()
        }
    }
}

impl<'a> AsRef<[u8]> for Ids<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}
