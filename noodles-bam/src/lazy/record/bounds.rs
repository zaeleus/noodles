use std::ops::{Range, RangeFrom};

pub const REFERENCE_SEQUENCE_ID_RANGE: Range<usize> = 0..4;
pub const ALIGNMENT_START_RANGE: Range<usize> = 4..8;
pub const MAPPING_QUALITY_RANGE: Range<usize> = 9..10;
pub const FLAGS_RANGE: Range<usize> = 14..16;
pub const MATE_REFERENCE_SEQUENCE_ID_RANGE: Range<usize> = 20..24;
pub const MATE_ALIGNMENT_START_RANGE: Range<usize> = 24..28;
pub const TEMPLATE_LENGTH_RANGE: Range<usize> = 28..32;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Bounds {
    pub read_name_end: usize,
    pub cigar_end: usize,
    pub sequence_end: usize,
    pub quality_scores_end: usize,
}

impl Bounds {
    pub fn read_name_range(&self) -> Range<usize> {
        TEMPLATE_LENGTH_RANGE.end..self.read_name_end
    }

    pub fn cigar_range(&self) -> Range<usize> {
        self.read_name_end..self.cigar_end
    }

    pub fn sequence_range(&self) -> Range<usize> {
        self.cigar_end..self.sequence_end
    }

    pub fn quality_scores_range(&self) -> Range<usize> {
        self.sequence_end..self.quality_scores_end
    }

    pub fn data_range(&self) -> RangeFrom<usize> {
        self.quality_scores_end..
    }
}
