mod iter;

use std::io;

use noodles_sam as sam;

use self::iter::Iter;
use super::Feature;

pub(super) struct Cigar<'r, 'c: 'r> {
    features: &'r [Feature<'c>],
    is_unmapped: bool,
    read_length: usize,
}

impl<'r, 'c: 'r> Cigar<'r, 'c> {
    pub(super) fn new(features: &'r [Feature<'c>], is_unmapped: bool, read_length: usize) -> Self {
        Self {
            features,
            is_unmapped,
            read_length,
        }
    }
}

impl sam::alignment::record::Cigar for Cigar<'_, '_> {
    fn is_empty(&self) -> bool {
        self.is_unmapped
    }

    fn len(&self) -> usize {
        if self.is_unmapped {
            0
        } else {
            self.iter().count()
        }
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<sam::alignment::record::cigar::Op>> + '_> {
        use std::iter;

        use sam::alignment::record::cigar::iter::TrySimplify;

        if self.is_unmapped {
            Box::new(iter::empty())
        } else {
            Box::new(TrySimplify::new(
                Iter::new(self.features, self.read_length).map(Ok),
            ))
        }
    }
}
