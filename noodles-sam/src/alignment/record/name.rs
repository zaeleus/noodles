/// An alignment record name.
pub trait Name {
    /// Returns the name as a byte slice.
    fn as_bytes(&self) -> &[u8];
}
