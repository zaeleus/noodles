use std::{
    io,
    ops::{Range, RangeFrom},
};

use bstr::ByteSlice;

#[derive(Clone, Eq, PartialEq)]
pub(super) struct Bounds {
    reference_sequence_name_end: usize,
    source_end: usize,
    type_end: usize,
    start_end: usize,
    end_end: usize,
    score_end: usize,
    strand_end: usize,
    phase_end: usize,
}

impl Bounds {
    pub(super) fn index(mut src: &[u8]) -> io::Result<Self> {
        let mut bounds = Self::default();
        let mut len = 0;

        len += read_required_field(&mut src)?;
        bounds.reference_sequence_name_end = len;

        len += read_required_field(&mut src)?;
        bounds.source_end = len;

        len += read_required_field(&mut src)?;
        bounds.type_end = len;

        len += read_required_field(&mut src)?;
        bounds.start_end = len;

        len += read_required_field(&mut src)?;
        bounds.end_end = len;

        len += read_required_field(&mut src)?;
        bounds.score_end = len;

        len += read_required_field(&mut src)?;
        bounds.strand_end = len;

        len += read_required_field(&mut src)?;
        bounds.phase_end = len;

        Ok(bounds)
    }

    pub fn reference_sequence_name_range(&self) -> Range<usize> {
        0..sans_delimiter(self.reference_sequence_name_end)
    }

    pub fn source_range(&self) -> Range<usize> {
        self.reference_sequence_name_end..sans_delimiter(self.source_end)
    }

    pub fn type_range(&self) -> Range<usize> {
        self.source_end..sans_delimiter(self.type_end)
    }

    pub fn start_range(&self) -> Range<usize> {
        self.type_end..sans_delimiter(self.start_end)
    }

    pub fn end_range(&self) -> Range<usize> {
        self.start_end..sans_delimiter(self.end_end)
    }

    pub fn score_range(&self) -> Range<usize> {
        self.end_end..sans_delimiter(self.score_end)
    }

    pub fn strand_range(&self) -> Range<usize> {
        self.score_end..sans_delimiter(self.strand_end)
    }

    pub fn phase_range(&self) -> Range<usize> {
        self.strand_end..sans_delimiter(self.phase_end)
    }

    pub fn attributes_range(&self) -> RangeFrom<usize> {
        self.phase_end..
    }
}

impl Default for Bounds {
    fn default() -> Self {
        Self {
            reference_sequence_name_end: 1,
            source_end: 2,
            type_end: 3,
            start_end: 4,
            end_end: 5,
            score_end: 6,
            strand_end: 7,
            phase_end: 8,
        }
    }
}

fn sans_delimiter(i: usize) -> usize {
    // SAFETY: Bound ends always include the delimiter.
    i - 1
}

fn read_required_field(src: &mut &[u8]) -> io::Result<usize> {
    let (len, is_eol) = read_field(src);

    if is_eol {
        Err(io::Error::from(io::ErrorKind::UnexpectedEof))
    } else {
        Ok(len)
    }
}

fn read_field(src: &mut &[u8]) -> (usize, bool) {
    const DELIMITER: u8 = b'\t';

    let (len, is_eol) = if let Some(i) = src.as_bstr().find_byte(DELIMITER) {
        (i + 1, false)
    } else {
        (src.len(), true)
    };

    *src = &src[len..];

    (len, is_eol)
}
