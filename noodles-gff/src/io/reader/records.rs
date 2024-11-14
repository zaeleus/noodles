use std::io::{self, BufRead};

use crate::{DirectiveBuf, LineBuf, Record};

use super::Lines;

/// Returns an iterator over records of a GFF reader.
///
/// This filters lines for only records. It stops at either EOF or when the `FASTA` directive is
/// read, whichever comes first.
///
/// This is created by calling [`crate::Reader::records`].
pub struct Records<'a, R> {
    lines: Lines<'a, R>,
}

impl<'a, R> Records<'a, R>
where
    R: BufRead,
{
    pub(crate) fn new(lines: Lines<'a, R>) -> Self {
        Self { lines }
    }
}

impl<'a, R> Iterator for Records<'a, R>
where
    R: BufRead,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.lines.next()? {
                Ok(line) => match line {
                    LineBuf::Directive(DirectiveBuf::StartOfFasta) => return None,
                    LineBuf::Record(r) => return Some(Ok(r)),
                    _ => {}
                },
                Err(e) => return Some(Err(e)),
            }
        }
    }
}
