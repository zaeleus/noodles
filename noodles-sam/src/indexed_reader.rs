//! Indexed SAM reader.

mod builder;

pub use self::builder::Builder;

use std::io::{self, Read, Seek};

use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi as csi;

use super::{alignment::Record, lazy, reader::Records, Header, Reader};

/// An indexed SAM reader.
pub struct IndexedReader<R> {
    inner: Reader<bgzf::Reader<R>>,
    index: csi::Index,
}

impl<R> IndexedReader<R>
where
    R: Read,
{
    /// Creates an indexed SAM reader.
    pub fn new(inner: R, index: csi::Index) -> Self {
        Self {
            inner: Reader::new(bgzf::Reader::new(inner)),
            index,
        }
    }

    /// Returns a reference to the underlying reader.
    pub fn get_ref(&self) -> &bgzf::Reader<R> {
        self.inner.get_ref()
    }

    /// Returns a mutable reference to the underlying reader.
    pub fn get_mut(&mut self) -> &mut bgzf::Reader<R> {
        self.inner.get_mut()
    }

    /// Returns the underlying reader.
    pub fn into_inner(self) -> bgzf::Reader<R> {
        self.inner.into_inner()
    }

    /// Reads the SAM header.
    pub fn read_header(&mut self) -> io::Result<Header> {
        self.inner.read_header()
    }

    /// Reads a single SAM record.
    pub fn read_record(&mut self, header: &Header, record: &mut Record) -> io::Result<usize> {
        self.inner.read_record(header, record)
    }

    /// Returns an iterator over records starting from the current stream position.
    pub fn records<'a>(&'a mut self, header: &'a Header) -> Records<'a, bgzf::Reader<R>> {
        self.inner.records(header)
    }

    /// Reads a single record without eagerly decoding its fields.
    pub fn read_lazy_record(&mut self, record: &mut lazy::Record) -> io::Result<usize> {
        self.inner.read_lazy_record(record)
    }

    /// Returns the associated index.
    pub fn index(&self) -> &csi::Index {
        &self.index
    }
}

impl<R> IndexedReader<R>
where
    R: Read + Seek,
{
    /// Returns an iterator over records that intersect the given region.
    ///
    /// To query for unmapped records, use [`Self::query_unmapped`].
    pub fn query<'a>(
        &'a mut self,
        header: &'a Header,
        region: &Region,
    ) -> io::Result<impl Iterator<Item = io::Result<Record>> + 'a> {
        self.inner.query(header, &self.index, region)
    }

    /// Returns an iterator of unmapped records after querying for the unmapped region.
    pub fn query_unmapped<'a>(
        &'a mut self,
        header: &'a Header,
    ) -> io::Result<impl Iterator<Item = io::Result<Record>> + 'a> {
        self.inner.query_unmapped(header, &self.index)
    }
}
