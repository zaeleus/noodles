use std::io::{self, BufRead};

use crate::{Record, record::Sequence};

use super::Reader;

/// An iterator over records of a FASTA reader.
///
/// This is created by calling [`Reader::records`].
pub struct Records<'a, R> {
    inner: &'a mut Reader<R>,
    line_buf: String,
    sequence_buf: Vec<u8>,
}

impl<'a, R> Records<'a, R>
where
    R: BufRead,
{
    pub(crate) fn new(inner: &'a mut Reader<R>) -> Self {
        Self {
            inner,
            line_buf: String::new(),
            sequence_buf: Vec::new(),
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

        self.sequence_buf.clear();

        match self.inner.read_sequence(&mut self.sequence_buf) {
            Ok(_) => {
                // Move the filled buffer into Sequence (zero-copy).
                // Replace with a new Vec pre-allocated to the same capacity
                // so the next record's read_to_end doesn't need to grow.
                let capacity = self.sequence_buf.capacity();
                let sequence = Sequence::from(std::mem::replace(
                    &mut self.sequence_buf,
                    Vec::with_capacity(capacity),
                ));
                let record = Record::new(definition, sequence);
                Some(Ok(record))
            }
            Err(e) => Some(Err(e)),
        }
    }
}
