use std::io;

use noodles_bgzf as bgzf;

use super::Query;
use crate::Record;

pub struct Records<'r, R> {
    inner: Query<'r, R>,
    record: Record,
}

impl<'r, R> Records<'r, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    pub(super) fn new(inner: Query<'r, R>) -> Self {
        Self {
            inner,
            record: Record::default(),
        }
    }
}

impl<R> Iterator for Records<'_, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.read_record(&mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
