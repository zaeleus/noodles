use std::io::{self, BufRead};

use crate::record::data::Field;

use super::Reader;

/// An iterator over BAM record data fields.
///
/// This is created by calling [`bam::record::data::Reader::fields`].
///
/// [`bam::record::data::Reader::fields`]: struct.Reader.html#method.fields
pub struct Fields<R>
where
    R: BufRead,
{
    reader: Reader<R>,
}

impl<R> Fields<R>
where
    R: BufRead,
{
    pub(crate) fn new(reader: Reader<R>) -> Self {
        Self { reader }
    }
}

impl<R> Iterator for Fields<R>
where
    R: BufRead,
{
    type Item = io::Result<Field>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_field() {
            Ok(Some(field)) => Some(Ok(field)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}
