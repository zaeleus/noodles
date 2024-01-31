use std::io;

/// Variant record info field array value values.
#[allow(clippy::len_without_is_empty)]
pub trait Values<'a, N> {
    /// Returns the number of values.
    fn len(&self) -> usize;

    /// Returns an iterator over values.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<N>>> + '_>;
}

const DELIMITER: char = ',';
const MISSING: char = '.';

impl<'a> Values<'a, char> for &'a str {
    fn len(&self) -> usize {
        (*self).len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<char>>> + '_> {
        Box::new(
            self.split(DELIMITER)
                .flat_map(|t| t.chars())
                .map(|c| match c {
                    MISSING => None,
                    _ => Some(c),
                })
                .map(Ok),
        )
    }
}
