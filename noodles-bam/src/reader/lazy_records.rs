use std::io::{self, Read};

use super::Reader;
use crate::lazy;

/// An iterator over lazily-evalulated records of a BAM reader.
///
/// This is created by calling [`Reader::lazy_records`].
pub struct LazyRecords<'a, R> {
    reader: &'a mut Reader<R>,
    record: lazy::Record,
}

impl<'a, R> LazyRecords<'a, R>
where
    R: Read,
{
    pub(super) fn new(reader: &'a mut Reader<R>) -> Self {
        Self {
            reader,
            record: lazy::Record::default(),
        }
    }
}

impl<'a, R> Iterator for LazyRecords<'a, R>
where
    R: Read,
{
    type Item = io::Result<lazy::Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_lazy_record(&mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
