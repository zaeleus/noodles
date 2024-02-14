/// Variant record filters.
pub trait Filters {
    /// Returns whether there are any filters.
    fn is_empty(&self) -> bool;

    /// Returns the number of filters.
    fn len(&self) -> usize;

    /// Returns an iterator over filters.
    fn iter(&self) -> Box<dyn Iterator<Item = &str> + '_>;

    /// Returns whether this is a `PASS` filter.
    fn is_pass(&self) -> bool {
        const PASS: &str = "PASS";

        let mut filters = self.iter();

        if let Some(filter) = filters.next() {
            filter == PASS && filters.next().is_none()
        } else {
            false
        }
    }
}
