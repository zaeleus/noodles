use std::ops::{Range, RangeFrom};

pub(super) const REFERENCE_SEQUENCE_ID_RANGE: Range<usize> = 0..4;
pub(super) const POSITION_RANGE: Range<usize> = 4..8;
pub(super) const SPAN_RANGE: Range<usize> = 8..12;
pub(super) const QUALITY_SCORE_RANGE: Range<usize> = 12..16;
pub(super) const INFO_FIELD_COUNT_RANGE: Range<usize> = 16..18;
pub(super) const ALLELE_COUNT_RANGE: Range<usize> = 18..20;
pub(super) const SAMPLE_COUNT_RANGE: Range<usize> = 20..23;
pub(super) const FORMAT_KEY_COUNT_INDEX: usize = 23;

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Bounds {
    pub(super) ids_range: Range<usize>,
    pub(super) reference_bases_range: Range<usize>,
    pub(super) alternate_bases_end: usize,
    pub(super) filters_end: usize,
}

impl Bounds {
    pub(super) fn ids_range(&self) -> Range<usize> {
        self.ids_range.clone()
    }

    pub(super) fn reference_bases_range(&self) -> Range<usize> {
        self.reference_bases_range.clone()
    }

    pub(super) fn alternate_bases_range(&self) -> Range<usize> {
        self.reference_bases_range.end..self.alternate_bases_end
    }

    pub(super) fn filters_range(&self) -> Range<usize> {
        self.alternate_bases_end..self.filters_end
    }

    pub(super) fn info_range(&self) -> RangeFrom<usize> {
        self.filters_end..
    }
}
