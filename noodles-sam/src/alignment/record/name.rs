/// An alignment record name.
pub trait Name {
    /// Returns the name as a byte slice.
    fn as_bytes(&self) -> &[u8];
}

impl Name for Box<dyn Name + '_> {
    fn as_bytes(&self) -> &[u8] {
        (**self).as_bytes()
    }
}
