use std::io::{self, BufRead};

use super::Reader;
use crate::{Header, Record};

/// An iterator over records of a VCF reader.
///
/// This is created by calling [`Reader::records`].
pub struct Records<'r, 'h, R> {
    inner: &'r mut Reader<R>,
    header: &'h Header,
    line_buf: String,
}

impl<'r, 'h, R> Records<'r, 'h, R>
where
    R: BufRead,
{
    pub(crate) fn new(inner: &'r mut Reader<R>, header: &'h Header) -> Self {
        Self {
            inner,
            header,
            line_buf: String::new(),
        }
    }
}

impl<'r, 'h, R> Iterator for Records<'r, 'h, R>
where
    R: BufRead,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        self.line_buf.clear();

        match self.inner.read_record(&mut self.line_buf) {
            Ok(0) => None,
            Ok(_) => Some(
                Record::try_from_str(&self.line_buf, self.header)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            ),
            Err(e) => Some(Err(e)),
        }
    }
}
