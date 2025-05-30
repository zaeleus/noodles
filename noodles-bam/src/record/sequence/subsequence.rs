use noodles_sam as sam;

use super::{Iter, decode_base};

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

    /// Returns the base at the given index.
    pub fn get(&self, i: usize) -> Option<u8> {
        let j = self.start + i;

        if j < self.end {
            let k = j / 2;
            let b = self.src[k];

            if j % 2 == 0 {
                Some(decode_base(b >> 4))
            } else {
                Some(decode_base(b))
            }
        } else {
            None
        }
    }

    /// Splits the subsequence into two subsequences at the given index.
    pub fn split_at_checked(&self, mid: usize) -> Option<(Self, Self)> {
        let mid = self.start + mid;

        if mid <= self.end {
            Some((
                Self::new(self.src, self.start, mid),
                Self::new(self.src, mid, self.end),
            ))
        } else {
            None
        }
    }

    /// Returns an iterator over the bases in the subsequence.
    pub fn iter(&self) -> impl Iterator<Item = u8> + '_ {
        Iter::new(self.src, self.start, self.end)
    }
}

impl sam::alignment::record::Sequence for Subsequence<'_> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn get(&self, i: usize) -> Option<u8> {
        self.get(i)
    }

    fn split_at_checked(
        &self,
        mid: usize,
    ) -> Option<(
        Box<dyn sam::alignment::record::Sequence + '_>,
        Box<dyn sam::alignment::record::Sequence + '_>,
    )> {
        self.split_at_checked(mid).map(|(left, right)| {
            (
                Box::new(left) as Box<dyn sam::alignment::record::Sequence + '_>,
                Box::new(right) as Box<dyn sam::alignment::record::Sequence + '_>,
            )
        })
    }

    fn iter(&self) -> Box<dyn Iterator<Item = u8> + '_> {
        Box::new(self.iter())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SRC: &[u8] = &[0x12];

    #[test]
    fn test_is_empty() {
        let subsequence = Subsequence::new(SRC, 0, 0);
        assert!(subsequence.is_empty());

        let subsequence = Subsequence::new(SRC, 0, 2);
        assert!(!subsequence.is_empty());
    }

    #[test]
    fn test_len() {
        let subsequence = Subsequence::new(SRC, 0, 0);
        assert_eq!(subsequence.len(), 0);

        let subsequence = Subsequence::new(SRC, 0, 2);
        assert_eq!(subsequence.len(), 2);
    }

    #[test]
    fn test_get() {
        let subsequence = Subsequence::new(SRC, 0, 0);
        assert!(subsequence.get(0).is_none());

        let subsequence = Subsequence::new(SRC, 1, 2);
        assert_eq!(subsequence.get(0), Some(b'C'));
    }

    #[test]
    fn test_split_at_checked() {
        let src = [0x12, 0x48];
        let subsequence = Subsequence::new(&src, 1, 3);

        assert_eq!(
            subsequence.split_at_checked(1),
            Some((Subsequence::new(&src, 1, 2), Subsequence::new(&src, 2, 3)))
        );

        assert!(subsequence.split_at_checked(4).is_none());
    }

    #[test]
    fn test_iter() {
        let subsequence = Subsequence::new(SRC, 0, 0);
        assert!(subsequence.iter().next().is_none());

        let subsequence = Subsequence::new(SRC, 0, 2);
        assert_eq!(subsequence.iter().collect::<Vec<_>>(), [b'A', b'C']);
    }
}
