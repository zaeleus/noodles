#![allow(dead_code)]

mod indexed_raw_reader;
mod raw_reader;

use std::io::{self, BufRead, Seek};

use noodles_core::Region;

use crate::Record;

pub use self::{indexed_raw_reader::IndexedRawReader, raw_reader::RawReader};

pub enum Inner<R> {
    Raw(RawReader<R>),
    IndexedRaw(IndexedRawReader<R>),
}

impl<R> Inner<R> {
    pub fn get_ref(&self) -> &R {
        match self {
            Self::Raw(inner) => inner.get_ref(),
            Self::IndexedRaw(inner) => inner.get_ref(),
        }
    }

    pub fn get_mut(&mut self) -> &mut R {
        match self {
            Self::Raw(inner) => inner.get_mut(),
            Self::IndexedRaw(inner) => inner.get_mut(),
        }
    }

    pub fn into_inner(self) -> R {
        match self {
            Self::Raw(inner) => inner.into_inner(),
            Self::IndexedRaw(inner) => inner.into_inner(),
        }
    }
}

impl<R> Inner<R>
where
    R: BufRead,
{
    pub fn read_definition(&mut self, buf: &mut String) -> io::Result<usize> {
        match self {
            Self::Raw(inner) => inner.read_definition(buf),
            Self::IndexedRaw(inner) => inner.read_definition(buf),
        }
    }

    pub fn read_sequence(&mut self, buf: &mut Vec<u8>) -> io::Result<usize> {
        match self {
            Self::Raw(inner) => inner.read_sequence(buf),
            Self::IndexedRaw(inner) => inner.read_sequence(buf),
        }
    }
}

impl<R> Inner<R>
where
    R: BufRead + Seek,
{
    pub fn query(&mut self, region: &Region) -> io::Result<Record> {
        match self {
            Self::Raw(_) => Err(io::Error::from(io::ErrorKind::Unsupported)),
            Self::IndexedRaw(inner) => inner.query(region),
        }
    }
}
