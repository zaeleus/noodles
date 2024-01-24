mod iter;

use noodles_sam as sam;

use self::iter::Iter;

/// A BAM record sequence.
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
        Iter::new(self.as_ref(), self.len())
    }
}

impl<'a> sam::alignment::record::Sequence for Sequence<'a> {
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
