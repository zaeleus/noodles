use std::io;

use crate::Header;

/// Variant record filters.
pub trait Filters {
    /// Returns whether there are any filters.
    fn is_empty(&self) -> bool;

    /// Returns the number of filters.
    fn len(&self) -> usize;

    /// Returns an iterator over filters.
    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<&'a str>> + 'a>;

    /// Returns whether this is a `PASS` filter.
    fn is_pass(&self, header: &Header) -> io::Result<bool> {
        const PASS: &str = "PASS";

        let mut filters = self.iter(header);

        if let Some(result) = filters.next() {
            let filter = result?;
            Ok(filter == PASS && filters.next().is_none())
        } else {
            Ok(false)
        }
    }
}
