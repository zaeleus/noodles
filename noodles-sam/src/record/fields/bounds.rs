use std::ops::{Range, RangeFrom};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Bounds {
    pub(crate) name_end: usize,
    pub(crate) flags_end: usize,
    pub(crate) reference_sequence_name_end: usize,
    pub(crate) alignment_start_end: usize,
    pub(crate) mapping_quality_end: usize,
    pub(crate) cigar_end: usize,
    pub(crate) mate_reference_sequence_name_end: usize,
    pub(crate) mate_alignment_start_end: usize,
    pub(crate) template_length_end: usize,
    pub(crate) sequence_end: usize,
    pub(crate) quality_scores_end: usize,
}

impl Bounds {
    pub fn name_range(&self) -> Range<usize> {
        0..self.name_end
    }

    pub fn flags_range(&self) -> Range<usize> {
        self.name_end..self.flags_end
    }

    pub fn reference_sequence_name_range(&self) -> Range<usize> {
        self.flags_end..self.reference_sequence_name_end
    }

    pub fn alignment_start_range(&self) -> Range<usize> {
        self.reference_sequence_name_end..self.alignment_start_end
    }

    pub fn mapping_quality_range(&self) -> Range<usize> {
        self.alignment_start_end..self.mapping_quality_end
    }

    pub fn cigar_range(&self) -> Range<usize> {
        self.mapping_quality_end..self.cigar_end
    }

    pub fn mate_reference_sequence_name_range(&self) -> Range<usize> {
        self.cigar_end..self.mate_reference_sequence_name_end
    }

    pub fn mate_alignment_start_range(&self) -> Range<usize> {
        self.mate_reference_sequence_name_end..self.mate_alignment_start_end
    }

    pub fn template_length_range(&self) -> Range<usize> {
        self.mate_alignment_start_end..self.template_length_end
    }

    pub fn sequence_range(&self) -> Range<usize> {
        self.template_length_end..self.sequence_end
    }

    pub fn quality_scores_range(&self) -> Range<usize> {
        self.sequence_end..self.quality_scores_end
    }

    pub fn data_range(&self) -> RangeFrom<usize> {
        self.quality_scores_end..
    }
}

impl Default for Bounds {
    fn default() -> Self {
        Self {
            name_end: 1,
            flags_end: 2,
            reference_sequence_name_end: 3,
            alignment_start_end: 4,
            mapping_quality_end: 7,
            cigar_end: 8,
            mate_reference_sequence_name_end: 9,
            mate_alignment_start_end: 10,
            template_length_end: 11,
            sequence_end: 12,
            quality_scores_end: 13,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ranges() {
        let bounds = Bounds::default();
        assert_eq!(bounds.name_range(), 0..1);
        assert_eq!(bounds.flags_range(), 1..2);
        assert_eq!(bounds.reference_sequence_name_range(), 2..3);
        assert_eq!(bounds.alignment_start_range(), 3..4);
        assert_eq!(bounds.mapping_quality_range(), 4..7);
        assert_eq!(bounds.cigar_range(), 7..8);
        assert_eq!(bounds.mate_reference_sequence_name_range(), 8..9);
        assert_eq!(bounds.mate_alignment_start_range(), 9..10);
        assert_eq!(bounds.template_length_range(), 10..11);
        assert_eq!(bounds.sequence_range(), 11..12);
        assert_eq!(bounds.quality_scores_range(), 12..13);
        assert_eq!(bounds.data_range(), 13..);
    }
}
