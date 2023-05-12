use std::io;

use noodles_core::Region;

use super::IndexedRecord;

/// An iterator that filters indexed records that intersect the given region.
pub struct FilterByRegion<'r, I, R>
where
    I: Iterator<Item = io::Result<R>>,
    R: IndexedRecord,
{
    records: I,
    region: &'r Region,
}

impl<'r, I, R> FilterByRegion<'r, I, R>
where
    I: Iterator<Item = io::Result<R>>,
    R: IndexedRecord,
{
    /// Creates a filtered indexed records iterator.
    pub fn new(records: I, region: &'r Region) -> Self {
        Self { records, region }
    }
}

impl<'r, I, R> Iterator for FilterByRegion<'r, I, R>
where
    I: Iterator<Item = io::Result<R>>,
    R: IndexedRecord,
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

fn intersects<R>(record: &R, region: &Region) -> bool
where
    R: IndexedRecord,
{
    record.indexed_reference_sequence_name() == region.name()
        && record.indexed_interval().intersects(region.interval())
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position;

    #[test]
    fn test_next() -> Result<(), Box<dyn std::error::Error>> {
        impl IndexedRecord for (&str, Position, Position) {
            fn indexed_reference_sequence_name(&self) -> &str {
                self.0
            }

            fn indexed_start_position(&self) -> Position {
                self.1
            }

            fn indexed_end_position(&self) -> Position {
                self.2
            }
        }

        let records = [
            Ok(("sq0", Position::try_from(5)?, Position::try_from(13)?)),
            Ok(("sq1", Position::try_from(5)?, Position::try_from(13)?)),
            Ok(("sq1", Position::try_from(21)?, Position::try_from(55)?)),
            Ok(("sq1", Position::try_from(34)?, Position::try_from(55)?)),
            Ok(("sq1", Position::try_from(89)?, Position::try_from(144)?)),
            Ok(("sq2", Position::try_from(5)?, Position::try_from(13)?)),
        ];

        let region = "sq1:21-55".parse()?;
        let actual: Vec<_> =
            FilterByRegion::new(records.into_iter(), &region).collect::<Result<_, _>>()?;

        let expected = [
            ("sq1", Position::try_from(21)?, Position::try_from(55)?),
            ("sq1", Position::try_from(34)?, Position::try_from(55)?),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
