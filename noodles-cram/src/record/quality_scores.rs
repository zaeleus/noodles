mod iter;

use std::io;

use noodles_sam as sam;

use self::iter::Iter;
use super::Feature;

pub struct QualityScores<'r> {
    features: &'r [Feature],
    read_length: usize,
}

impl<'r> QualityScores<'r> {
    pub fn new(features: &'r [Feature], read_length: usize) -> Self {
        Self {
            features,
            read_length,
        }
    }
}

impl sam::alignment::record::QualityScores for QualityScores<'_> {
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn len(&self) -> usize {
        self.read_length
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u8>> + '_> {
        Box::new(Iter::new(self.features, self.read_length))
    }
}
