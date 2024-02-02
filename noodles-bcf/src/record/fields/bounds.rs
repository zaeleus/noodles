use std::ops::Range;

pub(super) const REFERENCE_SEQUENCE_ID_RANGE: Range<usize> = 0..4;
pub(super) const POSITION_RANGE: Range<usize> = 4..8;
pub(super) const SPAN_RANGE: Range<usize> = 8..12;
pub(super) const QUALITY_SCORE_RANGE: Range<usize> = 12..16;
pub(super) const SAMPLE_COUNT_RANGE: Range<usize> = 20..23;
pub(super) const FORMAT_KEY_COUNT_INDEX: usize = 23;

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Bounds {
    pub(super) ids_range: Range<usize>,
    pub(super) reference_bases_range: Range<usize>,
}

impl Bounds {
    pub(super) fn ids_range(&self) -> Range<usize> {
        self.ids_range.clone()
    }

    pub(super) fn reference_bases_range(&self) -> Range<usize> {
        self.reference_bases_range.clone()
    }
}
