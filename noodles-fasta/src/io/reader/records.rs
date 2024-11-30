use std::io::{self, BufRead};

use crate::{record::Sequence, Record};

use super::Reader;

/// An iterator over records of a FASTA reader.
///
/// This is created by calling [`Reader::records`].
pub struct Records<'a, R> {
    inner: &'a mut Reader<R>,
    line_buf: String,
}

impl<'a, R> Records<'a, R>
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

impl<R> Iterator for Records<'_, R>
where
    R: BufRead,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        self.line_buf.clear();

        match self.inner.read_definition(&mut self.line_buf) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(e) => return Some(Err(e)),
        }

        let definition = match self.line_buf.parse() {
            Ok(d) => d,
            Err(e) => return Some(Err(io::Error::new(io::ErrorKind::InvalidData, e))),
        };

        let mut sequence_buf = Vec::new();

        match self.inner.read_sequence(&mut sequence_buf) {
            Ok(_) => {
                let record = Record::new(definition, Sequence::from(sequence_buf));
                Some(Ok(record))
            }
            Err(e) => Some(Err(e)),
        }
    }
}
