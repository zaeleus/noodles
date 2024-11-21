use std::io::{self, BufRead};

use crate::{directive_buf::key, DirectiveBuf, LineBuf, RecordBuf};

use super::LineBufs;

/// Returns an iterator over records of a GFF reader.
///
/// This filters lines for only records. It stops at either EOF or when the `FASTA` directive is
/// read, whichever comes first.
///
/// This is created by calling [`crate::Reader::records`].
pub struct RecordBufs<'a, R> {
    lines: LineBufs<'a, R>,
}

impl<'a, R> RecordBufs<'a, R>
where
    R: BufRead,
{
    pub(crate) fn new(lines: LineBufs<'a, R>) -> Self {
        Self { lines }
    }
}

impl<'a, R> Iterator for RecordBufs<'a, R>
where
    R: BufRead,
{
    type Item = io::Result<RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.lines.next()? {
                Ok(line) => match line {
                    LineBuf::Directive(DirectiveBuf::Other(key, _))
                        if key == key::START_OF_FASTA =>
                    {
                        return None
                    }
                    LineBuf::Record(r) => return Some(Ok(r)),
                    _ => {}
                },
                Err(e) => return Some(Err(e)),
            }
        }
    }
}
