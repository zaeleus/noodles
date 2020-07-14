use super::{Base, Sequence};

/// An iterator over bases of a sequence.
///
/// This is created by calling [`Sequence::bases`].
///
/// [`Sequence::bases`]: struct.Sequence.html#method.bases
pub struct Bases<'a> {
    sequence: &'a Sequence<'a>,
    head: usize,
    tail: usize,
    remaining: usize,
}

impl<'a> Bases<'a> {
    pub(crate) fn new(sequence: &'a Sequence<'_>) -> Self {
        let tail = if sequence.is_empty() {
            0
        } else {
            sequence.n_chars() - 1
        };

        Self {
            sequence,
            head: 0,
            tail,
            remaining: sequence.n_chars(),
        }
    }
}

impl<'a> Iterator for Bases<'a> {
    type Item = Base;

    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            return None;
        }

        let symbol = self.sequence[self.head];
        self.head += 1;
        self.remaining -= 1;
        Some(symbol)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}

impl<'a> DoubleEndedIterator for Bases<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            return None;
        }

        let symbol = self.sequence[self.tail];
        self.tail -= 1;
        self.remaining -= 1;
        Some(symbol)
    }
}
