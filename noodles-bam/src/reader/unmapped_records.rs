use std::io::{self, Read};

use crate::Record;

use super::Reader;

/// An iterator over unmapped records of a BAM reader.
///
/// This is created by calling [`bam::Reader::query_unmapped`].
///
/// [`bam::Reader::query_unmapped`]: struct.Reader.html#method.query_unmapped
pub struct UnmappedRecords<'a, R>
where
    R: Read,
{
    reader: &'a mut Reader<R>,
    record: Record,
}

impl<'a, R> UnmappedRecords<'a, R>
where
    R: Read,
{
    pub(crate) fn new(reader: &'a mut Reader<R>) -> Self {
        Self {
            reader,
            record: Record::default(),
        }
    }
}

impl<'a, R> Iterator for UnmappedRecords<'a, R>
where
    R: Read,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.reader.read_record(&mut self.record) {
                Ok(0) => return None,
                Ok(_) => {
                    if self.record.flags().is_unmapped() {
                        return Some(Ok(self.record.clone()));
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }
    }
}
