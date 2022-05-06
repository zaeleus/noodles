use std::io::{self, BufRead};

use crate::Record;

use super::{record::Fields, Reader};

/// An iterator over records of a SAM reader.
///
/// This is created by calling [`Reader::records`].
pub struct Records<'a, R> {
    inner: &'a mut Reader<R>,
    fields: Fields,
    line_buf: String,
}

impl<'a, R> Records<'a, R>
where
    R: BufRead,
{
    pub(crate) fn new(inner: &'a mut Reader<R>, fields: Fields) -> Self {
        Self {
            inner,
            fields,
            line_buf: String::new(),
        }
    }
}

impl<'a, R> Iterator for Records<'a, R>
where
    R: BufRead,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        self.line_buf.clear();

        match self.inner.read_record(&mut self.line_buf) {
            Ok(0) => None,
            Ok(_) => Some(
                Record::parse_with_fields(&self.line_buf, self.fields)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            ),
            Err(e) => Some(Err(e)),
        }
    }
}
