//! Indexed BAM reader.

mod builder;

use std::io::{self, Read};

use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;
use noodles_sam::{self as sam, alignment::RecordBuf};

pub use self::builder::Builder;
use super::{
    Reader,
    reader::{Query, RecordBufs, Records},
};
use crate::Record;

/// An indexed BAM reader.
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
    R: Read,
{
    /// Reads the SAM header.
    pub fn read_header(&mut self) -> io::Result<sam::Header> {
        self.inner.read_header()
    }

    /// Reads a record into an alignment record buffer.
    pub fn read_record_buf(
        &mut self,
        header: &sam::Header,
        record: &mut RecordBuf,
    ) -> io::Result<usize> {
        self.inner.read_record_buf(header, record)
    }

    /// Reads a record.
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        self.inner.read_record(record)
    }

    /// Returns an iterator over alignment record buffers starting from the current stream
    /// position.
    pub fn record_bufs<'r, 'h: 'r>(&'r mut self, header: &'h sam::Header) -> RecordBufs<'r, 'h, R> {
        self.inner.record_bufs(header)
    }

    /// Returns an iterator over records.
    pub fn records(&mut self) -> Records<'_, R> {
        self.inner.records()
    }

    /// Returns the associated index.
    pub fn index(&self) -> &dyn BinningIndex {
        &self.index
    }
}

impl<R> IndexedReader<bgzf::io::Reader<R>>
where
    R: Read,
{
    /// Creates an indexed BAM reader.
    pub fn new<I>(inner: R, index: I) -> Self
    where
        I: BinningIndex + 'static,
    {
        Self {
            inner: Reader::new(inner),
            index: Box::new(index),
        }
    }
}

impl<R> IndexedReader<R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    /// Returns an iterator over records that intersect the given region.
    pub fn query<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
        region: &Region,
    ) -> io::Result<Query<'r, R>> {
        self.inner.query(header, &self.index, region)
    }

    /// Returns an iterator of unmapped records after querying for the unmapped region.
    pub fn query_unmapped(&mut self) -> io::Result<impl Iterator<Item = io::Result<Record>>> {
        self.inner.query_unmapped(&self.index)
    }
}
