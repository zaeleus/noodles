use std::io::{self, BufRead, Seek};

use noodles_core::Region;

use super::raw_reader::RawReader;
use crate::{fai, Record};

pub struct IndexedRawReader<R> {
    inner: RawReader<R>,
    index: fai::Index,
}

impl<R> IndexedRawReader<R>
where
    R: BufRead,
{
    pub fn new(inner: R, index: fai::Index) -> Self {
        Self {
            inner: RawReader::new(inner),
            index,
        }
    }
}

impl<R> IndexedRawReader<R> {
    pub fn get_ref(&self) -> &R {
        self.inner.get_ref()
    }

    pub fn get_mut(&mut self) -> &mut R {
        self.inner.get_mut()
    }

    pub fn into_inner(self) -> R {
        self.inner.into_inner()
    }
}

impl<R> IndexedRawReader<R>
where
    R: BufRead,
{
    pub fn read_definition(&mut self, buf: &mut String) -> io::Result<usize> {
        self.inner.read_definition(buf)
    }

    pub fn read_sequence(&mut self, buf: &mut Vec<u8>) -> io::Result<usize> {
        self.inner.read_sequence(buf)
    }
}

impl<R> IndexedRawReader<R>
where
    R: BufRead + Seek,
{
    pub fn query(&mut self, region: &Region) -> io::Result<Record> {
        self.inner.query(&self.index, region)
    }
}
