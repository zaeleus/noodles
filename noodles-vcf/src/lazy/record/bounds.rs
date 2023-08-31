use std::ops::{Range, RangeFrom};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Bounds {
    pub chromosome_end: usize,
    pub position_end: usize,
    pub ids_end: usize,
    pub reference_bases_end: usize,
    pub alternate_bases_end: usize,
    pub quality_score_end: usize,
    pub filters_end: usize,
    pub info_end: usize,
}

impl Bounds {
    pub fn chromosome_range(&self) -> Range<usize> {
        0..self.chromosome_end
    }

    pub fn position_range(&self) -> Range<usize> {
        self.chromosome_end..self.position_end
    }

    pub fn ids_range(&self) -> Range<usize> {
        self.position_end..self.ids_end
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
        self.quality_score_end..self.filters_end
    }

    pub fn genotypes_range(&self) -> RangeFrom<usize> {
        self.info_end..
    }
}
