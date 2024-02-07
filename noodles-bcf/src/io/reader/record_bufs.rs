use std::io::{self, Read};

use noodles_vcf::{self as vcf, variant::RecordBuf};

use super::Reader;

/// An iterator over records of a BCF reader.
///
/// This is created by calling [`Reader::records`].
pub struct RecordBufs<'r, 'h, R>
where
    R: Read,
{
    reader: &'r mut Reader<R>,
    header: &'h vcf::Header,
    record: RecordBuf,
}

impl<'r, 'h, R> RecordBufs<'r, 'h, R>
where
    R: Read,
{
    pub(crate) fn new(reader: &'r mut Reader<R>, header: &'h vcf::Header) -> Self {
        Self {
            reader,
            header,
            record: RecordBuf::default(),
        }
    }
}

impl<'r, 'h, R> Iterator for RecordBufs<'r, 'h, R>
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
