use std::io::{self, BufRead};

use crate::Record;

use super::Reader;

/// An iterator over records of a FASTQ reader.
///
/// This is created by calling [`Reader::records`].
pub struct Records<'a, R> {
    inner: &'a mut Reader<R>,
    buf: Record,
}

impl<'a, R> Records<'a, R>
where
    R: BufRead,
{
    pub(crate) fn new(inner: &'a mut Reader<R>) -> Self {
        Self {
            inner,
            buf: Record::default(),
        }
    }
}

impl<'a, R> Iterator for Records<'a, R>
where
    R: BufRead,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        self.buf.clear();

        match self.inner.read_record(&mut self.buf) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.buf.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
