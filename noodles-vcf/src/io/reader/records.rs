use std::io::{self, BufRead};

use super::Reader;
use crate::{Header, Record};

/// An iterator over records of a VCF reader.
///
/// This is created by calling [`Reader::records`].
pub struct Records<'r, 'h, R> {
    inner: &'r mut Reader<R>,
    header: &'h Header,
    record: Record,
}

impl<'r, 'h, R> Records<'r, 'h, R>
where
    R: BufRead,
{
    pub(crate) fn new(inner: &'r mut Reader<R>, header: &'h Header) -> Self {
        Self {
            inner,
            header,
            record: Record::default(),
        }
    }
}

impl<'r, 'h, R> Iterator for Records<'r, 'h, R>
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
