use std::io::{self, BufRead};

use super::LineBufs;
use crate::{directive_buf::key, feature::RecordBuf, LineBuf};

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

impl<R> Iterator for RecordBufs<'_, R>
where
    R: BufRead,
{
    type Item = io::Result<RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.lines.next()? {
                Ok(line) => match line {
                    LineBuf::Directive(directive) if directive.key() == key::FASTA => {
                        return None;
                    }
                    LineBuf::Record(r) => return Some(Ok(r)),
                    _ => {}
                },
                Err(e) => return Some(Err(e)),
            }
        }
    }
}
