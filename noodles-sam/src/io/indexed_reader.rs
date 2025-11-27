//! Indexed SAM reader.

mod builder;

use std::io::{self, Read, Seek};

use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi::BinningIndex;

pub use self::builder::Builder;
use super::{Reader, reader::RecordBufs};
use crate::{Header, Record, alignment::RecordBuf};

/// An indexed SAM reader.
pub struct IndexedReader<R> {
    inner: Reader<bgzf::io::Reader<R>>,
    index: Box<dyn BinningIndex>,
}

impl<R> IndexedReader<R> {
    /// Returns a reference to the underlying reader.
    pub fn get_ref(&self) -> &bgzf::io::Reader<R> {
        self.inner.get_ref()
    }

    /// Returns a mutable reference to the underlying reader.
    pub fn get_mut(&mut self) -> &mut bgzf::io::Reader<R> {
        self.inner.get_mut()
    }

    /// Returns the underlying reader.
    pub fn into_inner(self) -> bgzf::io::Reader<R> {
        self.inner.into_inner()
    }
}

impl<R> IndexedReader<R>
where
    R: Read,
{
    /// Creates an indexed SAM reader.
    pub fn new<I>(inner: R, index: I) -> Self
    where
        I: BinningIndex + 'static,
    {
        Self {
            inner: Reader::new(bgzf::io::Reader::new(inner)),
            index: Box::new(index),
        }
    }

    /// Reads the SAM header.
    pub fn read_header(&mut self) -> io::Result<Header> {
        self.inner.read_header()
    }

    /// Reads a record into an alignment record buffer.
    pub fn read_record_buf(
        &mut self,
        header: &Header,
        record: &mut RecordBuf,
    ) -> io::Result<usize> {
        self.inner.read_record_buf(header, record)
    }

    /// Returns an iterator over alignment record buffers starting from the current stream
    /// position.
    pub fn record_bufs<'a>(
        &'a mut self,
        header: &'a Header,
    ) -> RecordBufs<'a, bgzf::io::Reader<R>> {
        self.inner.record_bufs(header)
    }

    /// Reads a record.
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
    R: Read + Seek,
{
    /// Returns an iterator over records that intersect the given region.
    pub fn query<'r, 'h: 'r>(
        &'r mut self,
        header: &'h Header,
        region: &Region,
    ) -> io::Result<impl Iterator<Item = io::Result<Record>> + use<'r, 'h, R>> {
        self.inner
            .query(header, &self.index, region)
            .map(|query| query.into_iter())
    }

    /// Returns an iterator of unmapped records after querying for the unmapped region.
    pub fn query_unmapped<'r>(
        &'r mut self,
    ) -> io::Result<impl Iterator<Item = io::Result<Record>> + use<'r, R>> {
        self.inner.query_unmapped(&self.index)
    }
}
