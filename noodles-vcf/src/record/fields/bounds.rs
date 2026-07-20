use std::ops::{Range, RangeFrom};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Bounds {
    pub reference_sequence_name_end: usize,
    pub variant_start_end: usize,
    pub ids_end: usize,
    pub reference_bases_end: usize,
    pub alternate_bases_end: usize,
    pub quality_score_end: usize,
    pub filters_end: usize,
    pub info_end: usize,
}

impl Bounds {
    pub fn reference_sequence_name_range(&self) -> Range<usize> {
        0..self.reference_sequence_name_end
    }

    pub fn variant_start_range(&self) -> Range<usize> {
        self.reference_sequence_name_end..self.variant_start_end
    }

    pub fn ids_range(&self) -> Range<usize> {
        self.variant_start_end..self.ids_end
    }

    pub fn reference_bases_range(&self) -> Range<usize> {
        self.ids_end..self.reference_bases_end
    }

    pub fn alternate_bases_range(&self) -> Range<usize> {
        self.reference_bases_end..self.alternate_bases_end
    }

    pub fn quality_score_range(&self) -> Range<usize> {
        self.alternate_bases_end..self.quality_score_end
    }

    pub fn filters_range(&self) -> Range<usize> {
        self.quality_score_end..self.filters_end
    }

    pub fn info_range(&self) -> Range<usize> {
        self.filters_end..self.info_end
    }

    pub fn samples_range(&self) -> RangeFrom<usize> {
        self.info_end..
    }
}

impl Default for Bounds {
    fn default() -> Self {
        Self {
            reference_sequence_name_end: 1,
            variant_start_end: 2,
            ids_end: 3,
            reference_bases_end: 4,
            alternate_bases_end: 5,
            quality_score_end: 6,
            filters_end: 7,
            info_end: 8,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ranges() {
        let bounds = Bounds::default();
        assert_eq!(bounds.reference_sequence_name_range(), 0..1);
        assert_eq!(bounds.variant_start_range(), 1..2);
        assert_eq!(bounds.ids_range(), 2..3);
        assert_eq!(bounds.reference_bases_range(), 3..4);
        assert_eq!(bounds.alternate_bases_range(), 4..5);
        assert_eq!(bounds.quality_score_range(), 5..6);
        assert_eq!(bounds.filters_range(), 6..7);
        assert_eq!(bounds.info_range(), 7..8);
        assert_eq!(bounds.samples_range(), 8..);
    }
}
