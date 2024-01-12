use std::io::{self, BufRead};

use super::Reader;
use crate::{alignment::RecordBuf, Header};

/// An iterator over record buffers of a SAM reader.
///
/// This is created by calling [`Reader::records`].
pub struct RecordBufs<'a, R> {
    inner: &'a mut Reader<R>,
    header: &'a Header,
    record: RecordBuf,
}

impl<'a, R> RecordBufs<'a, R>
where
    R: BufRead,
{
    pub(crate) fn new(inner: &'a mut Reader<R>, header: &'a Header) -> Self {
        Self {
            inner,
            header,
            record: RecordBuf::default(),
        }
    }
}

impl<'a, R> Iterator for RecordBufs<'a, R>
where
    R: BufRead,
{
    type Item = io::Result<RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.read_record_buf(self.header, &mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
