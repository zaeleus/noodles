use std::io::{self, BufRead};

use super::Reader;
use crate::{
    Record,
    record::{Definition, Sequence},
};

/// An iterator over records of a FASTA reader.
///
/// This is created by calling [`Reader::records`].
pub struct Records<'a, R> {
    inner: &'a mut Reader<R>,
    definition_buf: Definition,
}

impl<'a, R> Records<'a, R>
where
    R: BufRead,
{
    pub(crate) fn new(inner: &'a mut Reader<R>) -> Self {
        Self {
            inner,
            definition_buf: Definition::default(),
        }
    }
}

impl<R> Iterator for Records<'_, R>
where
    R: BufRead,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.read_definition(&mut self.definition_buf) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(e) => return Some(Err(e)),
        }

        let mut sequence_buf = Vec::new();

        match self.inner.read_sequence(&mut sequence_buf) {
            Ok(_) => {
                let record = Record::new(self.definition_buf.clone(), Sequence::from(sequence_buf));
                Some(Ok(record))
            }
            Err(e) => Some(Err(e)),
        }
    }
}
