use std::io;

use noodles_core::Region;

use super::Record;

/// A FASTA index.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Index(Vec<Record>);

impl Index {
    /// Returns start position of the given region.
    pub fn query(&self, region: &Region) -> io::Result<u64> {
        self.as_ref()
            .iter()
            .find(|record| record.name() == region.name())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid reference sequence name: {}", region.name()),
                )
            })
            .and_then(|record| record.query(region.interval()))
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
