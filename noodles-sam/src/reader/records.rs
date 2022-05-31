use std::io::{self, BufRead};

use super::Reader;
use crate::{alignment::Record, Header};

/// An iterator over records of a SAM reader.
///
/// This is created by calling [`Reader::records`].
pub struct Records<'a, R> {
    inner: &'a mut Reader<R>,
    header: &'a Header,
    record: Record,
}

impl<'a, R> Records<'a, R>
where
    R: BufRead,
{
    pub(crate) fn new(inner: &'a mut Reader<R>, header: &'a Header) -> Self {
        Self {
            inner,
            header,
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
        match self.inner.read_record(self.header, &mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
