use std::{array, iter::FusedIterator, slice};

use super::decode_base;

/// A BAM record sequence bases iterator.
pub struct Iter<'a> {
    iter: slice::Iter<'a, u8>,
    front: Option<array::IntoIter<u8, 2>>,
    back: Option<array::IntoIter<u8, 2>>,
}

impl<'a> Iter<'a> {
    // [start, end)
    pub(super) fn new(bases: &'a [u8], start: usize, end: usize) -> Self {
        let i = start / 2;
        let j = end.div_ceil(2);
        let mut iter = bases[i..j].iter();

        let front = if start % 2 == 0 {
            None
        } else {
            iter.next().map(|&n| discard_front_decoded_bases(n))
        };

        let base_count = end - start;

        // This assumes `bases.len() * 2` is only ever `base_count` or `base_count` + 1.
        let back = if bases.len() * 2 > base_count {
            iter.next_back().map(|&n| discard_back_decoded_bases(n))
        } else {
            None
        };

        Self { iter, front, back }
    }
}

impl Iterator for Iter<'_> {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(base) = self.front.as_mut().and_then(|iter| iter.next()) {
                return Some(base);
            }

            self.front = match self.iter.next() {
                Some(n) => Some(decoded_bases(*n)),
                None => return self.back.as_mut().and_then(|iter| iter.next()),
            };
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let mut n = self.iter.len() * 2;

        if let Some(ref iter) = self.front {
            n += iter.len();
        }

        if let Some(ref iter) = self.back {
            n += iter.len();
        }

        (n, Some(n))
    }
}

impl DoubleEndedIterator for Iter<'_> {
    fn next_back(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(base) = self.back.as_mut().and_then(|iter| iter.next_back()) {
                return Some(base);
            }

            self.back = match self.iter.next_back() {
                Some(n) => Some(decoded_bases(*n)),
                None => return self.front.as_mut().and_then(|iter| iter.next()),
            };
        }
    }
}

impl ExactSizeIterator for Iter<'_> {}

impl FusedIterator for Iter<'_> {}

pub(super) fn decoded_bases(n: u8) -> array::IntoIter<u8, 2> {
    [decode_base(n >> 4), decode_base(n)].into_iter()
}

fn discard_front_decoded_bases(n: u8) -> array::IntoIter<u8, 2> {
    let mut bases = decoded_bases(n);
    bases.next();
    bases
}

fn discard_back_decoded_bases(n: u8) -> array::IntoIter<u8, 2> {
    let mut bases = decoded_bases(n);
    bases.next_back();
    bases
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_next() {
        let mut iter = Iter::new(&[], 0, 0);
        assert!(iter.next().is_none());

        let iter = Iter::new(&[0x12, 0x40], 0, 3);
        assert_eq!(iter.collect::<Vec<_>>(), [b'A', b'C', b'G']);

        let iter = Iter::new(&[0x12, 0x48], 0, 4);
        assert_eq!(iter.collect::<Vec<_>>(), [b'A', b'C', b'G', b'T']);

        let iter = Iter::new(&[0x12, 0x48], 1, 3);
        assert_eq!(iter.collect::<Vec<_>>(), [b'C', b'G']);
    }

    #[test]
    fn test_size_hint() {
        fn t(mut iter: Iter, mut len: usize) {
            assert_eq!(iter.size_hint(), (len, Some(len)));

            while iter.next().is_some() {
                len -= 1;
                assert_eq!(iter.size_hint(), (len, Some(len)));
            }
        }

        t(Iter::new(&[], 0, 0), 0);
        t(Iter::new(&[0x12, 0x48], 0, 3), 3);
        t(Iter::new(&[0x12, 0x48], 0, 4), 4);
        t(Iter::new(&[0x12, 0x48], 1, 3), 2);
    }

    #[test]
    fn test_next_back() {
        let mut iter = Iter::new(&[], 0, 0);
        assert!(iter.next_back().is_none());

        let iter = Iter::new(&[0x12, 0x40], 0, 3);
        assert_eq!(iter.rev().collect::<Vec<_>>(), [b'G', b'C', b'A']);

        let iter = Iter::new(&[0x12, 0x48], 0, 4);
        assert_eq!(iter.rev().collect::<Vec<_>>(), [b'T', b'G', b'C', b'A']);

        let iter = Iter::new(&[0x12, 0x48], 1, 3);
        assert_eq!(iter.rev().collect::<Vec<_>>(), [b'G', b'C']);
    }
}
