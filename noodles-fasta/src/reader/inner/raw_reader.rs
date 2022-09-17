use std::io::{self, BufRead, Seek, SeekFrom};

use noodles_core::Region;

use crate::{fai, Record};

pub struct RawReader<R> {
    inner: R,
}

impl<R> RawReader<R>
where
    R: BufRead,
{
    pub fn new(inner: R) -> Self {
        Self { inner }
    }
}

impl<R> RawReader<R> {
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> RawReader<R>
where
    R: BufRead,
{
    pub fn read_definition(&mut self, buf: &mut String) -> io::Result<usize> {
        crate::reader::read_line(self.get_mut(), buf)
    }

    pub fn read_sequence(&mut self, buf: &mut Vec<u8>) -> io::Result<usize> {
        crate::reader::read_sequence(self.get_mut(), buf)
    }
}

impl<R> RawReader<R>
where
    R: BufRead + Seek,
{
    pub fn query(&mut self, index: &[fai::Record], region: &Region) -> io::Result<Record> {
        use crate::{
            reader::interval_to_slice_range,
            record::{Definition, Sequence},
        };

        let index_record = index
            .iter()
            .find(|record| record.name() == region.name())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid reference sequence name: {}", region.name()),
                )
            })?;

        let pos = index_record.offset();
        self.get_mut().seek(SeekFrom::Start(pos))?;

        let definition = Definition::new(region.to_string(), None);

        let mut raw_sequence = Vec::new();
        self.read_sequence(&mut raw_sequence)?;

        let range = interval_to_slice_range(region.interval(), raw_sequence.len());
        let sequence = Sequence::from(raw_sequence[range].to_vec());

        Ok(Record::new(definition, sequence))
    }
}
