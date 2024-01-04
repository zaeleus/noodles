use std::io::{self, Read};

use super::Reader;
use crate::Record;

/// An iterator over lazily-evalulated records of a BAM reader.
///
/// This is created by calling [`Reader::lazy_records`].
pub struct LazyRecords<'a, R> {
    reader: &'a mut Reader<R>,
    record: Record,
}

impl<'a, R> LazyRecords<'a, R>
where
    R: Read,
{
    pub(super) fn new(reader: &'a mut Reader<R>) -> Self {
        Self {
            reader,
            record: Record::default(),
        }
    }
}

impl<'a, R> Iterator for LazyRecords<'a, R>
where
    R: Read,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_lazy_record(&mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
