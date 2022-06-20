use std::io;

use noodles_sam as sam;

/// Raw BAM record CIGAR operations.
#[derive(Debug, Eq, PartialEq)]
pub struct Cigar<'a>(&'a [u8]);

impl<'a> Cigar<'a> {
    pub(super) fn new(src: &'a [u8]) -> Self {
        Self(src)
    }

    /// Returns whether there are any CIGAR operations.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of CIGAR operations.
    ///
    /// This is _not_ the length of the buffer.
    pub fn len(&self) -> usize {
        self.0.len() / 4
    }
}

impl<'a> AsRef<[u8]> for Cigar<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> TryFrom<Cigar<'a>> for sam::record::Cigar {
    type Error = io::Error;

    fn try_from(bam_cigar: Cigar<'a>) -> Result<Self, Self::Error> {
        use crate::reader::record::get_cigar;

        let mut src = bam_cigar.0;
        let mut cigar = Self::default();
        let op_count = bam_cigar.len();
        get_cigar(&mut src, &mut cigar, op_count)?;

        Ok(cigar)
    }
}
