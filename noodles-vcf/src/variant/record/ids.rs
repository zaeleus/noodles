/// Variant record IDs.
pub trait Ids {
    /// Returns whethere there are any IDs.
    fn is_empty(&self) -> bool;

    /// Returns the number of IDs.
    fn len(&self) -> usize;

    /// Returns an iterator over IDs.
    fn iter(&self) -> Box<dyn Iterator<Item = &str> + '_>;
}
