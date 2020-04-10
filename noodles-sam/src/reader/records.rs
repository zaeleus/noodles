use std::io::{self, BufRead};

use crate::Record;

use super::Reader;

pub struct Records<'a, R> {
    inner: &'a mut Reader<R>,
    line_buf: String,
}

impl<'a, R> Records<'a, R>
where
    R: BufRead,
{
    pub fn new(inner: &'a mut Reader<R>) -> Self {
        Self {
            inner,
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
                self.line_buf
                    .parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            ),
            Err(e) => Some(Err(e)),
        }
    }
}
