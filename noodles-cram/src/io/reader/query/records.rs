use std::io::{self, Read, Seek};

use noodles_sam as sam;

use super::Query;

pub struct Records<'r, 'h: 'r, 'i: 'r, R>
where
    R: Read + Seek,
{
    inner: Query<'r, 'h, 'i, R>,
    record: sam::alignment::RecordBuf,
}

impl<'r, 'h: 'r, 'i: 'r, R> Records<'r, 'h, 'i, R>
where
    R: Read + Seek,
{
    pub(super) fn new(inner: Query<'r, 'h, 'i, R>) -> Self {
        Self {
            inner,
            record: sam::alignment::RecordBuf::default(),
        }
    }
}

impl<'r, 'h: 'r, 'i: 'r, R> Iterator for Records<'r, 'h, 'i, R>
where
    R: Read + Seek,
{
    type Item = io::Result<sam::alignment::RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.read_record_buf(&mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
