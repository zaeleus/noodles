use std::ops::{Range, RangeFrom};

#[derive(Default)]
pub(crate) struct Bounds {
    pub(crate) reference_sequence_name_end: usize,
    pub(crate) source_end: usize,
    pub(crate) type_end: usize,
    pub(crate) start_end: usize,
    pub(crate) end_end: usize,
    pub(crate) score_end: usize,
    pub(crate) strand_end: usize,
    pub(crate) phase_end: usize,
}

impl Bounds {
    pub fn reference_sequence_name_range(&self) -> Range<usize> {
        0..self.reference_sequence_name_end
    }

    pub fn source_range(&self) -> Range<usize> {
        self.reference_sequence_name_end..self.source_end
    }

    pub fn type_range(&self) -> Range<usize> {
        self.source_end..self.type_end
    }

    pub fn start_range(&self) -> Range<usize> {
        self.type_end..self.start_end
    }

    pub fn end_range(&self) -> Range<usize> {
        self.start_end..self.end_end
    }

    pub fn score_range(&self) -> Range<usize> {
        self.end_end..self.score_end
    }

    pub fn strand_range(&self) -> Range<usize> {
        self.score_end..self.strand_end
    }

    pub fn phase_range(&self) -> Range<usize> {
        self.strand_end..self.phase_end
    }

    pub fn attributes_range(&self) -> RangeFrom<usize> {
        self.phase_end..
    }
}
