use std::io::{self, BufRead};

use crate::LineBuf;

use super::Reader;

/// An iterator over lines of a GFF reader.
///
/// When using this, the caller is responsible to stop reading at either EOF or when the `FASTA`
/// directive is read, whichever comes first.
///
/// This is created by calling [`Reader::lines`].
pub struct LineBufs<'a, R> {
    inner: &'a mut Reader<R>,
    line_buf: String,
}

impl<'a, R> LineBufs<'a, R>
where
    R: BufRead,
{
    pub(crate) fn new(inner: &'a mut Reader<R>) -> Self {
        Self {
            inner,
            line_buf: String::new(),
        }
    }
}

impl<R> Iterator for LineBufs<'_, R>
where
    R: BufRead,
{
    type Item = io::Result<LineBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        self.line_buf.clear();

        match self.inner.read_line_buf(&mut self.line_buf) {
            Ok(0) => None,
            Ok(_) => match self.line_buf.parse() {
                Ok(line) => Some(Ok(line)),
                Err(e) => Some(Err(io::Error::new(io::ErrorKind::InvalidData, e))),
            },
            Err(e) => Some(Err(e)),
        }
    }
}
