use std::io;

use noodles_core::Region;

use super::indexed_records::Record;

/// An iterator that filters raw records that intersect the given region.
pub struct FilterByRegion<'r, I> {
    records: I,
    region: &'r Region,
}

impl<'r, I> FilterByRegion<'r, I>
where
    I: Iterator<Item = io::Result<Record>>,
{
    /// Creates a filtered indexed records iterator.
    pub fn new(records: I, region: &'r Region) -> Self {
        Self { records, region }
    }
}

impl<'r, I> Iterator for FilterByRegion<'r, I>
where
    I: Iterator<Item = io::Result<Record>>,
{
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let record = match self.records.next()? {
                Ok(r) => r,
                Err(e) => return Some(Err(e)),
            };

            if intersects(&record, self.region) {
                return Some(Ok(record));
            }
        }
    }
}

fn intersects(record: &Record, region: &Region) -> bool {
    record.reference_sequence_name() == region.name()
        && record.interval().intersects(region.interval())
}
