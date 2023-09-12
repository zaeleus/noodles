/// Raw SAM record data.
#[derive(Debug, Eq, PartialEq)]
pub struct Data<'a>(&'a [u8]);

impl<'a> Data<'a> {
    pub(super) fn new(buf: &'a [u8]) -> Self {
        Self(buf)
    }

    /// Returns whether there are any data fields.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of data fields.
    pub fn len(&self) -> usize {
        self.0.len()
    }
}

impl<'a> AsRef<[u8]> for Data<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}
