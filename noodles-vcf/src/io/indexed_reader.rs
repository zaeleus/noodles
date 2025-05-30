//! Indexed VCF reader.

mod builder;

pub use self::builder::Builder;

use std::io::{self, BufRead, Read};

use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;

use super::{
    Reader,
    reader::{Query, RecordBufs},
};
use crate::{Header, Record, variant::RecordBuf};

/// An indexed VCF reader.
pub struct IndexedReader<R> {
    inner: Reader<R>,
    index: Box<dyn BinningIndex>,
}

impl<R> IndexedReader<R> {
    /// Returns a reference to the underlying reader.
    pub fn get_ref(&self) -> &R {
        self.inner.get_ref()
    }

    /// Returns a mutable reference to the underlying reader.
    pub fn get_mut(&mut self) -> &mut R {
        self.inner.get_mut()
    }

    /// Returns the underlying reader.
    pub fn into_inner(self) -> R {
        self.inner.into_inner()
    }
}

impl<R> IndexedReader<R>
where
    R: BufRead,
{
    /// Reads the VCF header.
    pub fn read_header(&mut self) -> io::Result<Header> {
        self.inner.read_header()
    }

    /// Reads a single raw VCF record.
    pub fn read_record_buf(
        &mut self,
        header: &Header,
        record: &mut RecordBuf,
    ) -> io::Result<usize> {
        self.inner.read_record_buf(header, record)
    }

    /// Returns an iterator over records starting from the current stream position.
    pub fn record_bufs<'r, 'h: 'r>(&'r mut self, header: &'h Header) -> RecordBufs<'r, 'h, R> {
        self.inner.record_bufs(header)
    }

    /// Reads a single record without eagerly parsing its fields.
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        self.inner.read_record(record)
    }

    /// Returns an iterator over records.
    pub fn records(&mut self) -> impl Iterator<Item = io::Result<Record>> {
        self.inner.records()
    }

    /// Returns the associated index.
    pub fn index(&self) -> &dyn BinningIndex {
        &self.index
    }
}

impl<R> IndexedReader<R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    /// Returns an iterator over records that intersects the given region.
    pub fn query<'r, 'h>(
        &'r mut self,
        header: &'h Header,
        region: &Region,
    ) -> io::Result<Query<'r, 'h, R>> {
        self.inner.query(header, &self.index, region)
    }
}

impl<R> IndexedReader<bgzf::io::Reader<R>>
where
    R: Read,
{
    /// Creates an indexed VCF reader.
    pub fn new<I>(inner: R, index: I) -> Self
    where
        I: BinningIndex + 'static,
    {
        Self {
            inner: Reader::new(bgzf::io::Reader::new(inner)),
            index: Box::new(index),
        }
    }
}
