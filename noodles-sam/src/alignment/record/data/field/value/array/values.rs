use std::io;

/// An alignment record data field array value values.
pub trait Values<'a, N> {
    /// Returns an iterator over values.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<N>> + '_>;
}
