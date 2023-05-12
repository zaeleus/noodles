mod record;

pub use self::record::Record;

use std::io::{self, BufRead, Lines};

use self::record::parse_record;
use crate::index::{header::format::CoordinateSystem, Header};

/// An iterator over indexed records.
pub struct IndexedRecords<R> {
    lines: Lines<R>,
    line_comment_prefix: char,
    reference_sequence_name_index: usize,
    start_position_index: usize,
    end_position_index: Option<usize>,
    coordinate_system: CoordinateSystem,
}

impl<R> IndexedRecords<R>
where
    R: BufRead,
{
    /// Creates an indexed records iterator.
    pub fn new(reader: R, header: &Header) -> Self {
        Self {
            lines: reader.lines(),
            line_comment_prefix: char::from(header.line_comment_prefix()),
            reference_sequence_name_index: header.reference_sequence_name_index(),
            start_position_index: header.start_position_index(),
            end_position_index: header.end_position_index(),
            coordinate_system: header.format().coordinate_system(),
        }
    }
}

impl<R> Iterator for IndexedRecords<R>
where
    R: BufRead,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line = match self.lines.next()? {
                Ok(s) => s,
                Err(e) => return Some(Err(e)),
            };

            if line.starts_with(self.line_comment_prefix) {
                continue;
            }

            let result = parse_record(
                line,
                self.reference_sequence_name_index,
                self.start_position_index,
                self.end_position_index,
                self.coordinate_system,
            )
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e));

            return Some(result);
        }
    }
}
