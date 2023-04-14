use std::io::{self, Read};

use noodles_vcf as vcf;

use crate::header::StringMaps;

use super::Reader;

/// An iterator over records of a BCF reader.
///
/// This is created by calling [`Reader::records`].
pub struct Records<'a, R>
where
    R: Read,
{
    reader: &'a mut Reader<R>,
    header: &'a vcf::Header,
    string_maps: &'a StringMaps,
    record: vcf::Record,
}

impl<'a, R> Records<'a, R>
where
    R: Read,
{
    pub(crate) fn new(
        reader: &'a mut Reader<R>,
        header: &'a vcf::Header,
        string_maps: &'a StringMaps,
    ) -> Records<'a, R> {
        Self {
            reader,
            header,
            string_maps,
            record: vcf::Record::default(),
        }
    }
}

impl<'a, R> Iterator for Records<'a, R>
where
    R: Read,
{
    type Item = io::Result<vcf::Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self
            .reader
            .read_record(self.header, self.string_maps, &mut self.record)
        {
            Ok(0) => None,
            Ok(_) => Some(Ok(self.record.clone())),
            Err(e) => Some(Err(e)),
        }
    }
}
