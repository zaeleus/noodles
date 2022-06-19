use std::io;

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
}

impl<'a> AsRef<[u8]> for Sequence<'a> {
    fn as_ref(&self) -> &[u8] {
        self.src
    }
}

impl<'a> TryFrom<Sequence<'a>> for sam::record::Sequence {
    type Error = io::Error;

    fn try_from(bam_sequence: Sequence<'a>) -> Result<Self, Self::Error> {
        use crate::reader::record::get_sequence;

        let mut src = bam_sequence.src;
        let mut sam_sequence = Self::default();
        let base_count = bam_sequence.base_count;
        get_sequence(&mut src, &mut sam_sequence, base_count)?;

        Ok(sam_sequence)
    }
}
