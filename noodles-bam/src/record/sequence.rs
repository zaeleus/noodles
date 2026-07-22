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
    pub(crate) fn new(src: &'a [u8], base_count: usize) -> Self {
        Self { src, base_count }
    }

    /// Returns a byte slice of encoded bases.
    pub fn as_bytes(&self) -> &'a [u8] {
        self.src
    }

    /// Returns whether there are any bases.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
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
            let n = self.src[j];
            let [l, r] = decode_bases(n);
            let b = if i.is_multiple_of(2) { l } else { r };
            Some(b)
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
        sequence.iter().collect()
    }
}

const CODES: [[u8; 2]; 256] = build_codes();

const fn build_codes() -> [[u8; 2]; 256] {
    const BASES: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

    let mut table = [[0u8; 2]; 256];
    let mut i = 0;

    while i < 256 {
        table[i] = [BASES[i >> 4], BASES[i & 0xf]];
        i += 1;
    }

    table
}

fn decode_bases(n: u8) -> [u8; 2] {
    CODES[usize::from(n)]
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
    fn test_from_sequence_for_sam_alignment_record_buf_sequence() {
        use noodles_sam::alignment::record_buf::Sequence as SequenceBuf;

        let sequence = Sequence::new(&[0x12, 0x40], 3);
        let actual = SequenceBuf::from(sequence);
        let expected = SequenceBuf::from(b"ACG");
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_build_codes() {
        let expected = *b"=ACMGRSVTWYHKDBN";

        for (n, base) in expected.into_iter().enumerate() {
            assert_eq!(CODES[n << 4][0], base);
            assert_eq!(CODES[n][1], base);
        }
    }
}
