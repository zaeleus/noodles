use std::io;

use noodles_bgzf as bgzf;

use super::Reader;
use crate::Record;

pub struct Iter<'r, R> {
    reader: Reader<'r, R>,
    record: Record,
}

impl<'r, R> Iter<'r, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    pub(super) fn new(reader: Reader<'r, R>) -> Self {
        Self {
            reader,
            record: Record::default(),
        }
    }
}

impl<R> Iterator for Iter<'_, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
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
