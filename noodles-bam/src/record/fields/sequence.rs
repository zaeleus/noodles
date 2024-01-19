use noodles_sam as sam;

/// A raw BAM record sequence.
#[derive(Debug, Eq, PartialEq)]
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

    /// Returns an iterator over the bases in the sequence.
    pub fn iter(&self) -> impl Iterator<Item = u8> + '_ {
        self.src
            .iter()
            .flat_map(|&b| [decode_base(b >> 4), decode_base(b)])
            .take(self.base_count)
    }
}

impl<'a> sam::alignment::record::field::Sequence for Sequence<'a> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = u8> + '_> {
        Box::new(self.iter())
    }
}

impl<'a> AsRef<[u8]> for Sequence<'a> {
    fn as_ref(&self) -> &[u8] {
        self.src
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
    fn test_iter() {
        let sequence = Sequence::new(&[], 0);
        assert!(sequence.iter().next().is_none());

        let sequence = Sequence::new(&[0x12, 0x40], 3);
        let actual: Vec<_> = sequence.iter().collect();
        assert_eq!(actual, b"ACG");

        let sequence = Sequence::new(&[0x12, 0x48], 4);
        let actual: Vec<_> = sequence.iter().collect();
        assert_eq!(actual, b"ACGT");
    }

    #[test]
    fn test_sam_alignment_record_sequence_iter() {
        fn t(src: &[u8], base_count: usize, expected: &[u8]) {
            let sequence: &dyn sam::alignment::record::field::Sequence =
                &Sequence::new(src, base_count);
            let actual: Vec<_> = sequence.iter().collect();
            assert_eq!(actual, expected);
        }

        t(&[], 0, &[]);
        t(&[0x12, 0x40], 3, b"ACG");
        t(&[0x12, 0x48], 4, b"ACGT");
    }
}
