/// A feature record attributes field array value.
pub trait Array {
    /// Returns an iterator over values.
    fn iter(&self) -> Box<dyn Iterator<Item = &str> + '_>;
}
