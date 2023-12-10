/// An alignment record read name.
pub trait ReadName {
    /// Returns the read name as a byte slice.
    fn as_bytes(&self) -> &[u8];
}
