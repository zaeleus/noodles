use std::ops::Range;

pub(super) const REFERENCE_SEQUENCE_ID_RANGE: Range<usize> = 0..4;
pub(super) const POSITION_RANGE: Range<usize> = 4..8;
pub(super) const SPAN_RANGE: Range<usize> = 8..12;
pub(super) const QUALITY_SCORE_RANGE: Range<usize> = 12..16;

#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct Bounds;
