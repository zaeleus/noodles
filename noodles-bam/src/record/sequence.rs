mod iter;
mod subsequence;

use std::fmt::{self, Write};

use noodles_sam as sam;

use self::iter::Iter;
pub use self::subsequence::Subsequence;

/// A BAM record sequence.
#[derive(Eq, PartialEq)]
pub struct Sequence<'a> {
    src: &'a [u8],
    base_count: usize,
}

impl<'a> Sequence<'a> {
    pub(super) fn new(src: &'a [u8], base_count: usize) -> Self {
        Self { src, base_count }
    }

    /// Returns whether there are any bases.
    pub fn is_empty(&self) -> bool {
        self.src.is_empty()
    }

    /// Returns the number of bases in the sequence.
    ///
    /// This is _not_ the length of the buffer.
    pub fn len(&self) -> usize {
        self.base_count
    }

    /// Returns the base at the given index.
    pub fn get(&self, i: usize) -> Option<u8> {
        if i < self.len() {
            let j = i / 2;
            let b = self.src[j];

            if i % 2 == 0 {
                Some(decode_base(b >> 4))
            } else {
                Some(decode_base(b))
            }
        } else {
            None
        }
    }

    /// Splits the sequence into two subsequences at the given index.
    ///
    /// The left split contains bases from `[0, mid)`; and the right, `[mid, len)`. If `mid` is
    /// > `len`, this returns `None`.
    pub fn split_at_checked(&'a self, mid: usize) -> Option<(Subsequence<'a>, Subsequence<'a>)> {
        if mid <= self.len() {
            Some((
                Subsequence::new(self.as_ref(), 0, mid),
                Subsequence::new(self.as_ref(), mid, self.len()),
            ))
        } else {
            None
        }
    }

    /// Returns an iterator over the bases in the sequence.
    pub fn iter(&self) -> Iter<'a> {
        Iter::new(self.src, 0, self.len())
    }
}

impl sam::alignment::record::Sequence for Sequence<'_> {
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

impl AsRef<[u8]> for Sequence<'_> {
    fn as_ref(&self) -> &[u8] {
        self.src
    }
}

impl fmt::Debug for Sequence<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        struct BasesFormat<'a>(&'a Sequence<'a>);

        impl fmt::Debug for BasesFormat<'_> {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_char('"')?;

                let BasesFormat(sequence) = self;

                for b in sequence.iter() {
                    f.write_char(char::from(b))?;
                }

                f.write_char('"')?;

                Ok(())
            }
        }

        f.debug_tuple("Sequence").field(&BasesFormat(self)).finish()
    }
}

impl<'a> From<Sequence<'a>> for sam::alignment::record_buf::Sequence {
    fn from(sequence: Sequence<'a>) -> Self {
        Self::from(sequence.as_ref().to_vec())
    }
}

fn decode_base(n: u8) -> u8 {
    match n & 0x0f {
        0 => b'=',
        1 => b'A',
        2 => b'C',
        3 => b'M',
        4 => b'G',
        5 => b'R',
        6 => b'S',
        7 => b'V',
        8 => b'T',
        9 => b'W',
        10 => b'Y',
        11 => b'H',
        12 => b'K',
        13 => b'D',
        14 => b'B',
        15 => b'N',
        _ => unreachable!(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        let sequence = Sequence::new(&[], 0);
        assert!(sequence.is_empty());

        let sequence = Sequence::new(&[0x12, 0x40], 3);
        assert!(!sequence.is_empty());
    }

    #[test]
    fn test_len() {
        let sequence = Sequence::new(&[], 0);
        assert_eq!(sequence.len(), 0);

        let sequence = Sequence::new(&[0x12, 0x40], 3);
        assert_eq!(sequence.len(), 3);
    }

    #[test]
    fn test_get() {
        let sequence = Sequence::new(&[0x12, 0x40], 3);
        assert_eq!(sequence.get(0), Some(b'A'));
        assert_eq!(sequence.get(1), Some(b'C'));
        assert_eq!(sequence.get(2), Some(b'G'));
        assert!(sequence.get(3).is_none());
    }

    #[test]
    fn test_split_at_checked() {
        let src = [0x10];
        let sequence = Sequence::new(&src, 1);

        assert_eq!(
            sequence.split_at_checked(0),
            Some((Subsequence::new(&src, 0, 0), Subsequence::new(&src, 0, 1)))
        );

        assert_eq!(
            sequence.split_at_checked(1),
            Some((Subsequence::new(&src, 0, 1), Subsequence::new(&src, 1, 1)))
        );

        assert!(sequence.split_at_checked(2).is_none());
    }

    #[test]
    fn test_debug_fmt() {
        let sequence = Sequence::new(&[], 0);
        assert_eq!(format!("{sequence:?}"), r#"Sequence("")"#);

        let sequence = Sequence::new(&[0x12, 0x40], 3);
        assert_eq!(format!("{sequence:?}"), r#"Sequence("ACG")"#);

        let sequence = Sequence::new(&[0x12, 0x48], 4);
        assert_eq!(format!("{sequence:?}"), r#"Sequence("ACGT")"#);
    }

    #[test]
    fn test_decode_base() {
        assert_eq!(decode_base(0), b'=');
        assert_eq!(decode_base(1), b'A');
        assert_eq!(decode_base(2), b'C');
        assert_eq!(decode_base(3), b'M');
        assert_eq!(decode_base(4), b'G');
        assert_eq!(decode_base(5), b'R');
        assert_eq!(decode_base(6), b'S');
        assert_eq!(decode_base(7), b'V');
        assert_eq!(decode_base(8), b'T');
        assert_eq!(decode_base(9), b'W');
        assert_eq!(decode_base(10), b'Y');
        assert_eq!(decode_base(11), b'H');
        assert_eq!(decode_base(12), b'K');
        assert_eq!(decode_base(13), b'D');
        assert_eq!(decode_base(14), b'B');
        assert_eq!(decode_base(15), b'N');
    }
}
