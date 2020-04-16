use std::io::{self, Read};

use crate::Record;

use super::Reader;

pub struct Records<'a, R: Read> {
    reader: &'a mut Reader<R>,
    record: Record,
}

impl<'a, R: Read> Records<'a, R> {
    pub fn new(reader: &'a mut Reader<R>) -> Records<R> {
        Self {
            reader,
            record: Record::default(),
        }
    }
}

impl<'a, R: Read> Iterator for Records<'a, R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record(&mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
