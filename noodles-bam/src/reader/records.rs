use std::io::{self, Read, Seek};

use crate::Record;

use super::Reader;

pub struct Records<'a, R: Read + Seek> {
    reader: &'a mut Reader<R>,
    record: Record,
}

impl<'a, R: Read + Seek> Records<'a, R> {
    pub fn new(reader: &'a mut Reader<R>) -> Records<R> {
        Self {
            reader,
            record: Record::new(),
        }
    }
}

impl<'a, R: Read + Seek> Iterator for Records<'a, R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record(&mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
