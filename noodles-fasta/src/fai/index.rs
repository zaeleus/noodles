use std::io;

use noodles_core::Region;

use super::Record;

/// A FASTA index.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Index(Vec<Record>);

impl Index {
    /// Returns start position of the given region.
    pub fn query(&self, region: &Region) -> io::Result<u64> {
        let record = self
            .as_ref()
            .iter()
            .find(|record| record.name() == region.name())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid reference sequence name: {}", region.name(),),
                )
            })?;

        let start = region
            .interval()
            .start()
            .map(|position| usize::from(position) - 1)
            .unwrap_or_default();

        let start =
            u64::try_from(start).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let pos = record.offset()
            + start / record.line_bases() * record.line_width()
            + start % record.line_bases();

        Ok(pos)
    }
}

impl AsRef<[Record]> for Index {
    fn as_ref(&self) -> &[Record] {
        &self.0
    }
}

impl From<Vec<Record>> for Index {
    fn from(records: Vec<Record>) -> Self {
        Self(records)
    }
}

impl From<Index> for Vec<Record> {
    fn from(index: Index) -> Self {
        index.0
    }
}
