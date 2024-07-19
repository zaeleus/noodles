use super::Iter;

/// A BAM record subsequence.
#[derive(Debug, Eq, PartialEq)]
pub struct Subsequence<'a> {
    src: &'a [u8],
    // [start, end)
    start: usize,
    end: usize,
}

impl<'a> Subsequence<'a> {
    pub(super) fn new(src: &'a [u8], start: usize, end: usize) -> Self {
        Self { src, start, end }
    }

    /// Returns whether there are any bases.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the number of bases.
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    /// Returns an iterator over the bases in the subsequence.
    pub fn iter(&self) -> impl Iterator<Item = u8> + '_ {
        Iter::new(self.src, self.start, self.end)
    }
}
