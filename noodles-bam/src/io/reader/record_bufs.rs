use std::io::{self, Read};

use noodles_sam::{self as sam, alignment::RecordBuf};

use super::Reader;

/// An iterator over records of a BAM reader.
///
/// This is created by calling [`Reader::record_bufs`].
pub struct RecordBufs<'r, 'h: 'r, R>
where
    R: Read,
{
    reader: &'r mut Reader<R>,
    header: &'h sam::Header,
    record: RecordBuf,
}

impl<'r, 'h: 'r, R> RecordBufs<'r, 'h, R>
where
    R: Read,
{
    pub(super) fn new(reader: &'r mut Reader<R>, header: &'h sam::Header) -> Self {
        Self {
            reader,
            header,
            record: RecordBuf::default(),
        }
    }
}

impl<R> Iterator for RecordBufs<'_, '_, R>
where
    R: Read,
{
    type Item = io::Result<RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record_buf(self.header, &mut self.record) {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
