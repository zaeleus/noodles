use std::io::{self, Read};

use crate::Record;

use super::Reader;

pub struct Records<'a, R>
where
    R: Read,
{
    reader: &'a mut Reader<R>,
    record: Record,
}

impl<'a, R> Records<'a, R>
where
    R: Read,
{
    pub fn new(reader: &'a mut Reader<R>) -> Records<'_, R> {
        Self {
            reader,
            record: Record::default(),
        }
    }
}

impl<'a, R> Iterator for Records<'a, R>
where
    R: Read,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record(&mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
