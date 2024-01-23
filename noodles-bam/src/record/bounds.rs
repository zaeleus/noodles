use std::ops::{Range, RangeFrom};

pub const REFERENCE_SEQUENCE_ID_RANGE: Range<usize> = 0..4;
pub const ALIGNMENT_START_RANGE: Range<usize> = 4..8;
pub const NAME_LENGTH_INDEX: usize = 8;
pub const MAPPING_QUALITY_INDEX: usize = 9;
pub const CIGAR_OP_COUNT_RANGE: Range<usize> = 12..14;
pub const FLAGS_RANGE: Range<usize> = 14..16;
pub const READ_LENGTH_RANGE: Range<usize> = 16..20;
pub const MATE_REFERENCE_SEQUENCE_ID_RANGE: Range<usize> = 20..24;
pub const MATE_ALIGNMENT_START_RANGE: Range<usize> = 24..28;
pub const TEMPLATE_LENGTH_RANGE: Range<usize> = 28..32;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Bounds {
    pub name_end: usize,
    pub cigar_end: usize,
    pub sequence_end: usize,
    pub quality_scores_end: usize,
}

impl Bounds {
    pub fn name_range(&self) -> Range<usize> {
        TEMPLATE_LENGTH_RANGE.end..self.name_end
    }

    pub fn cigar_range(&self) -> Range<usize> {
        self.name_end..self.cigar_end
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ranges() {
        let bounds = Bounds {
            name_end: 34,
            cigar_end: 38,
            sequence_end: 40,
            quality_scores_end: 42,
        };

        assert_eq!(bounds.name_range(), 32..34);
        assert_eq!(bounds.cigar_range(), 34..38);
        assert_eq!(bounds.sequence_range(), 38..40);
        assert_eq!(bounds.quality_scores_range(), 40..42);
        assert_eq!(bounds.data_range(), 42..);
    }
}
