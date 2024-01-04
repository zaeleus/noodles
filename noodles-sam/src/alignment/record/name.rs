/// An alignment record name.
#[allow(clippy::len_without_is_empty)]
pub trait Name {
    /// Returns the name as a byte slice.
    fn as_bytes(&self) -> &[u8];

    /// Returns the length.
    fn len(&self) -> usize {
        self.as_bytes().len()
    }
}

impl Name for Box<dyn Name + '_> {
    fn as_bytes(&self) -> &[u8] {
        (**self).as_bytes()
    }
}
