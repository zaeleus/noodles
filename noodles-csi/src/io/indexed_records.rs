mod record;

use noodles_core::Region;

pub use self::record::Record;

use std::io::{self, BufRead, Lines};

use self::record::parse_record;
use crate::index::{header::format::CoordinateSystem, Header};

use super::FilterByRegion;

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

    /// Creates an iterator that filters indexed records that intersect the given region.
    pub fn filter_by_region(self, region: &Region) -> FilterByRegion<Self, Record> {
        FilterByRegion::new(self, region)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_next() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"sq0\t8\t13
# noodles
sq0\t21\t34
";

        let reader = &data[..];

        let header = Header::builder()
            .set_start_position_index(1)
            .set_end_position_index(Some(2))
            .set_reference_sequence_names([String::from("sq0")].into_iter().collect())
            .build();

        let records: Vec<_> = IndexedRecords::new(reader, &header).collect::<Result<_, _>>()?;
        let lines: Vec<_> = records.iter().map(|r| r.as_ref()).collect();

        let expected = ["sq0\t8\t13", "sq0\t21\t34"];
        assert_eq!(lines, expected);

        Ok(())
    }
}
