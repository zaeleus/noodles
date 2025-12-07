use std::io;

use noodles_bgzf as bgzf;

use super::Query;
use crate::Record;

pub struct Records<'r, 'h: 'r, R> {
    inner: Query<'r, 'h, R>,
    record: Record,
}

impl<'r, 'h: 'r, R> Records<'r, 'h, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    pub(super) fn new(inner: Query<'r, 'h, R>) -> Self {
        Self {
            inner,
            record: Record::default(),
        }
    }
}

impl<R> Iterator for Records<'_, '_, R>
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
