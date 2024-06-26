use super::Record;

/// A FASTA index.
#[derive(Debug, Default, Eq, PartialEq)]
pub struct Index(Vec<Record>);

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
