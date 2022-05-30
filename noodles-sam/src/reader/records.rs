use std::io::{self, BufRead};

use super::Reader;
use crate::Record;

/// An iterator over records of a SAM reader.
///
/// This is created by calling [`Reader::records`].
pub struct Records<'a, R> {
    inner: &'a mut Reader<R>,
    record: Record,
}

impl<'a, R> Records<'a, R>
where
    R: BufRead,
{
    pub(crate) fn new(inner: &'a mut Reader<R>) -> Self {
        Self {
            inner,
            record: Record::default(),
        }
    }
}

impl<'a, R> Iterator for Records<'a, R>
where
    R: BufRead,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.read_record(&mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
