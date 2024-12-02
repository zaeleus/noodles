use std::io::{self, BufRead};

use super::Reader;
use crate::Line;

/// An iterator over lines of a GFF reader.
///
/// This is created by calling [`Reader::lines`].
pub struct Lines<'r, R> {
    reader: &'r mut Reader<R>,
    line: Line,
}

impl<'r, R> Lines<'r, R>
where
    R: BufRead,
{
    pub(super) fn new(reader: &'r mut Reader<R>) -> Self {
        Self {
            reader,
            line: Line::default(),
        }
    }
}

impl<R> Iterator for Lines<'_, R>
where
    R: BufRead,
{
    type Item = io::Result<Line>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_line(&mut self.line) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.line.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
